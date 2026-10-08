"""Recover locally available structures by exact names and descriptor checks.

This default workflow makes no network requests and leaves source files intact.
It writes a private local index for automatic association in the app.
"""
import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path
from bbb.chemistry import molecule_data
from bbb.library import read_public_libraries, source_name
from bbb.pipeline import canonical_row, read_file
from bbb.structures import STRUCTURE_FIELDS, identity_key

ROOT = Path(__file__).resolve().parent


def checks(row, structure):
    d = structure["descriptors"]
    comparisons = {}
    for field in STRUCTURE_FIELDS:
        try:
            original = float(row[field])
        except (TypeError, ValueError):
            return None
        tolerance = .02 if field == "MW" else .05 if field == "TPSA" else 0
        comparisons[field] = {"input": original, "rdkit": d[field],
                              "matches": abs(original-d[field]) <= tolerance+1e-9}
    return comparisons


def local_paper_candidates():
    path = ROOT/"references"/"jm9b01220_si_003.csv"
    candidates = defaultdict(list)
    if not path.exists():
        return candidates
    with path.open(newline="", encoding="utf-8-sig") as f:
        for line, values in enumerate(csv.reader(f), 1):
            if len(values) < 3 or not values[1].strip() or not values[2].strip():
                continue
            try:
                structure = molecule_data(values[1], False)
            except ValueError:
                continue
            candidates[values[2].strip().casefold()].append({
                "smiles": values[1], "structure": structure,
                "source_title": "Gupta 2019 supplementary structure data",
                "source_url": "https://doi.org/10.1021/acs.jmedchem.9b01220.s003",
                "source_locator": f"jm9b01220_si_003.csv line {line}",
                "match_method": "exact_name_and_descriptor_checks",
                "source_priority": 1,
            })
    return candidates


def recover(files, id_namespace='unknown'):
    candidates = local_paper_candidates()
    public = read_public_libraries()
    for name, values in public['by_name'].items():
        candidates[name].extend(values)
    records = {}
    unresolved = []
    input_rows = 0
    for path in files:
        for row_number, raw in enumerate(read_file(path.name, path.read_bytes()), 1):
            row = canonical_row(raw)
            input_rows += 1
            label = str(row["compound_id"] or "").strip()
            matches = []
            details = []
            if label.isdecimal():
                options = public['by_id'].get((id_namespace, label), []) if id_namespace != 'unknown' else []
                # Unknown namespaces stay unresolved. Candidate CIDs are exported
                # for review, never attached solely because digits coincide.
                review_options = public['by_id'].get(('pubchem', label), []) if id_namespace == 'unknown' else []
            else:
                # Bare paper labels such as 3a or s1 are not global identities.
                options = [] if re.fullmatch(r'\d+[a-z]{1,2}|[a-z]\d+', label, re.I) else candidates.get(source_name(label), [])
                review_options = []
            for candidate in [*options, *review_options]:
                try:
                    structure = candidate.get('structure') or molecule_data(candidate['smiles'], False)
                except ValueError as error:
                    details.append({'source': candidate['source_locator'], 'smiles': candidate['smiles'], 'error': str(error)})
                    continue
                candidate = {**candidate, 'structure': structure}
                comparison = checks(row, structure)
                details.append({'source': candidate['source_locator'], 'smiles': candidate['smiles'],
                                'source_compound_name': candidate.get('source_compound_name'), 'checks': comparison,
                                'identifier_namespace_confirmed': not label.isdecimal() or id_namespace != 'unknown'})
                # Exact identity from the named source, corroborated by MW, HA,
                # aromatic rings and TPSA. HBA definitions may differ by toolkit.
                if (not label.isdecimal() or id_namespace != 'unknown') and comparison and all(comparison[k]["matches"] for k in ("MW", "HA", "AroR", "TPSA", "HBD")):
                    matches.append((candidate, comparison))
            if matches:
                priority = min(c.get('source_priority', 1) for c, _ in matches)
                matches = [(c, comparison) for c, comparison in matches if c.get('source_priority', 1) == priority]
            graphs = {c["structure"]["canonical_smiles"] for c, _ in matches}
            key = identity_key(row)
            if key and len(graphs) == 1:
                c, comparison = matches[0]
                records[key] = {
                    "compound_id": label,
                    "input_descriptors": {k: row[k] for k in STRUCTURE_FIELDS},
                    "smiles": c["smiles"], "canonical_smiles": c["structure"]["canonical_smiles"],
                    "formula": c["structure"]["formula"],
                    "source_title": c["source_title"], "source_url": c["source_url"],
                    "source_locator": c["source_locator"],
                    "match_method": ('exact_' + id_namespace + '_id_and_descriptor_checks') if label.isdecimal() else c.get('match_method', 'exact_source_name_or_synonym_and_descriptor_checks'),
                    "checks": comparison, "allow_name_lookup": not label.isdecimal(),
                }
            else:
                reason = ('identifier_namespace_unconfirmed' if label.isdecimal() and id_namespace == 'unknown'
                          else 'ambiguous_source_structures' if len(graphs) > 1
                          else 'source_structure_or_descriptor_mismatch' if details
                          else 'absent_from_local_references')
                unresolved.append({"source_file": path.name, "row_id": row_number, "compound_id": label,
                                   'reason_code': reason,
                                   "reason": "No verified local structure match. Provide a source identifier or structure.",
                                   'candidates_for_review': details})
    return {"version": 1, "records": records}, {"input_rows": input_rows,
             "mapped_signatures": len(records), "unresolved_rows": len(unresolved),
             'identifier_namespace': id_namespace, 'public_libraries': public['datasets'], "unresolved": unresolved}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Match structures from local reference files")
    parser.add_argument("files", nargs="+", type=Path)
    parser.add_argument('--id-namespace', choices=['unknown', 'pubchem', 'gtopdb', 'chembl', 'drugcentral'], default='unknown',
                        help='Only set when the numeric ID namespace has been confirmed; default unknown')
    args = parser.parse_args()
    index, report = recover(args.files, args.id_namespace)
    output = ROOT/".local"
    output.mkdir(exist_ok=True)
    temporary = output/"structure_index.tmp.json"
    temporary.write_text(json.dumps(index, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    temporary.replace(output/"structure_index.json")
    (output/"structure_recovery_report.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in report.items() if k != "unresolved"}, ensure_ascii=False))
