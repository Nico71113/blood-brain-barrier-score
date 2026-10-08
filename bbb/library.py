"""Read public bulk files locally. This module never makes network requests."""
import csv
import html
import re
from collections import defaultdict
from pathlib import Path

REFERENCES = Path(__file__).resolve().parents[1] / 'references'


def source_name(value):
    # Presentation markup and trademark symbols carry no molecular identity.
    # Preserve punctuation and stereo prefixes; no fuzzy name matching.
    return html.unescape(re.sub(r'<[^>]*>', '', str(value))).replace('®', '').replace('™', '').strip().casefold()


def read_public_libraries(directory=REFERENCES):
    by_name, by_id, datasets = defaultdict(list), defaultdict(list), []
    path = directory / 'gtopdb_ligands.csv'
    if path.exists():
        with path.open(newline='', encoding='utf-8-sig') as stream:
            version = next(stream).strip().strip('"')
            rows = list(csv.DictReader(stream))
        datasets.append({'file': path.name, 'records': len(rows), 'version': version})
        for line, row in enumerate(rows, 3):
            if not row.get('SMILES'):
                continue
            candidate = {'smiles': row['SMILES'], 'source_priority': 2,
                         'source_title': 'IUPHAR/BPS Guide to PHARMACOLOGY',
                         'source_url': 'https://www.guidetopharmacology.org/DATA/ligands.csv',
                         'source_locator': f"gtopdb_ligands.csv line {line}; GtoPdb {row['Ligand ID']}",
                         'source_compound_name': row['Name']}
            names = {row['Name'], *row.get('Synonyms', '').split('|'), *row.get('INN', '').split('|')}
            for name in {source_name(n) for n in names} - {''}:
                by_name[name].append(candidate)
            for namespace, column in [('pubchem', 'PubChem CID'), ('gtopdb', 'Ligand ID'), ('chembl', 'ChEMBL ID')]:
                for identifier in row.get(column, '').split('|'):
                    if identifier.strip():
                        by_id[(namespace, identifier.strip())].append(candidate)
    path = directory / 'drugcentral_structures.smiles.tsv'
    if path.exists():
        with path.open(newline='', encoding='utf-8-sig') as stream:
            rows = list(csv.DictReader(stream, delimiter='\t'))
        datasets.append({'file': path.name, 'records': len(rows)})
        for line, row in enumerate(rows, 2):
            if not row.get('SMILES') or not row.get('INN'):
                continue
            candidate = {'smiles': row['SMILES'], 'source_priority': 3,
                         'source_title': 'DrugCentral',
                         'source_url': 'https://unmtid-dbs.net/download/DrugCentral/2021_09_01/structures.smiles.tsv',
                         'source_locator': f"drugcentral_structures.smiles.tsv line {line}; DrugCentral {row['ID']}",
                         'source_compound_name': row['INN']}
            by_name[source_name(row['INN'])].append(candidate)
            by_id[('drugcentral', row['ID'])].append(candidate)
    return {'by_name': by_name, 'by_id': by_id, 'datasets': datasets}
