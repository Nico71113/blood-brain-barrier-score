import csv
import io
import unittest
from unittest.mock import patch
from bbb.chemistry import molecule_data
from bbb.pipeline import analyze, batch_analyze, results_csv
from bbb.structures import identity_key, lookup_structure
from recover_structures import checks


class AssociationTests(unittest.TestCase):
    def setUp(self):
        self.molecule = molecule_data('CCO', False)
        self.row = dict(compound_id='local-id', **self.molecule['descriptors'], pKa=8)
        self.entry = dict(compound_id='local-id', smiles='CCO', canonical_smiles='CCO',
                          source_title='Local source fixture', source_url='https://example.invalid/source',
                          source_locator='row 1', match_method='exact_name_and_descriptor_checks',
                          checks=checks(self.row, self.molecule), allow_name_lookup=True)
        self.index = dict(version=1, records={identity_key(self.row): self.entry})

    def test_repeated_label_does_not_match_changed_descriptors(self):
        self.assertIs(lookup_structure(self.row, self.index), self.entry)
        self.assertIsNone(lookup_structure({**self.row, 'MW': 100}, self.index))
        self.assertIsNone(lookup_structure(dict(compound_id='local-id', MW=46.069), self.index))
        self.assertIsNone(identity_key({**self.row, 'HA': True}))
        self.assertEqual(identity_key(self.row), identity_key({**self.row, 'MW': str(self.row['MW'])}))

    def test_name_only_requires_one_verified_graph(self):
        self.assertIs(lookup_structure(dict(compound_id='LOCAL-ID'), self.index), self.entry)
        self.index['records']['other'] = {**self.entry, 'canonical_smiles': 'COC', 'smiles': 'COC'}
        self.assertIsNone(lookup_structure(dict(compound_id='local-id'), self.index))
        self.assertIsNone(lookup_structure(dict(compound_id='unknown'), self.index))

    def test_association_preserves_original_score_and_inputs(self):
        original = {**self.row, 'BBB SCORE': 3.2}
        with patch('bbb.pipeline.lookup_structure', side_effect=lambda r: lookup_structure(r, self.index)):
            result = analyze(original)
            summary = batch_analyze([original, {**self.row, 'MW': 100}])['summary']
        self.assertEqual(result['structure_origin'], 'local_association')
        self.assertEqual(result['structure']['formula'], 'C2H6O')
        self.assertEqual(result['bbb_score_input'], 3.2)
        self.assertEqual(result['descriptor_source'], 'uploaded_or_entered')
        self.assertIsNone(result['input']['smiles'])
        self.assertEqual(summary['recovered_structures'], 1)
        exported = next(csv.DictReader(io.StringIO(results_csv([result]))))
        self.assertEqual(exported['input_smiles'], 'CCO')
        self.assertEqual(exported['structure_source_locator'], 'row 1')

    def test_supplied_smiles_takes_priority(self):
        with patch('bbb.pipeline.lookup_structure') as lookup:
            result = analyze({**self.row, 'smiles': 'COC'})
            lookup.assert_not_called()
        self.assertEqual(result['structure_origin'], 'input_smiles')
        self.assertEqual(result['structure']['canonical_smiles'], 'COC')

    def test_source_metric_warning_only_for_original_metric_scoring(self):
        self.entry['checks']['HBA']['matches'] = False
        with patch('bbb.pipeline.lookup_structure', return_value=self.entry):
            by_name = analyze(dict(compound_id='local-id'))
            by_table = analyze(self.row)
        self.assertIsNone(by_name['bbb_score'])
        self.assertEqual(by_name['warnings'], [])
        self.assertTrue(any('HBA' in warning for warning in by_table['warnings']))


if __name__ == '__main__':
    unittest.main()
