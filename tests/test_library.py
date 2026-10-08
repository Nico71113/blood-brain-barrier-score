import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from bbb.chemistry import molecule_data
from bbb.library import source_name, read_public_libraries
from bbb.structures import identity_key
from recover_structures import recover


class LibraryTests(unittest.TestCase):
    def test_presentation_markup_preserves_stereo_and_punctuation(self):
        self.assertEqual(source_name('<i>Arimidex</i>&reg;'), 'arimidex')
        self.assertNotEqual(source_name('(S)-drug'), source_name('(R)-drug'))
        self.assertNotEqual(source_name('drug-1'), source_name('drug-2'))

    def test_public_aliases_and_namespaced_ids(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'gtopdb_ligands.csv'
            path.write_text('"# Test public dataset"\nLigand ID,Name,INN,Synonyms,SMILES,PubChem CID,ChEMBL ID\n9,ethanol,ethanol,trade&reg;|alcohol,CCO,702,CHEMBL545\n')
            library = read_public_libraries(Path(folder))
        self.assertEqual(len(library['by_name']['ethanol']), 1)
        self.assertEqual(library['by_name']['trade'][0]['smiles'], 'CCO')
        self.assertIn(('pubchem', '702'), library['by_id'])
        self.assertNotIn(('drugcentral', '702'), library['by_id'])

    def test_unconfirmed_numeric_ids_are_not_attached(self):
        row = dict(compound_id='702', **molecule_data('CCO', False)['descriptors'], pKa=8)
        candidate = dict(smiles='CCO', source_title='Public fixture', source_url='https://example.invalid/bulk',
                         source_locator='line 3', source_priority=2)
        library = dict(by_name={}, by_id={('pubchem', '702'): [candidate]}, datasets=[])
        with patch('recover_structures.local_paper_candidates', return_value={}), \
             patch('recover_structures.read_public_libraries', return_value=library), \
             patch('recover_structures.read_file', return_value=[row]):
            with tempfile.TemporaryDirectory() as folder:
                file = Path(folder) / 'table.csv'; file.write_bytes(b'fixture')
                unknown, report = recover([file])
                confirmed, _ = recover([file], 'pubchem')
                mismatch, _ = recover([file], 'gtopdb')
        self.assertEqual(unknown['records'], {})
        self.assertEqual(report['unresolved'][0]['reason_code'], 'identifier_namespace_unconfirmed')
        self.assertEqual(len(report['unresolved'][0]['candidates_for_review']), 1)
        self.assertIn(identity_key(row), confirmed['records'])
        self.assertFalse(confirmed['records'][identity_key(row)]['allow_name_lookup'])
        self.assertEqual(mismatch['records'], {})


if __name__ == '__main__':
    unittest.main()
