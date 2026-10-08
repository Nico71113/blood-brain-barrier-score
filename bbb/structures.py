"""Associate real, sourced structures with input records; no graph inference.

The local index is built by recover_structures.py. Compound labels alone are
not globally unique; a row association includes all six original descriptors.
"""
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

INDEX_PATH = Path(__file__).resolve().parents[1] / ".local" / "structure_index.json"
STRUCTURE_FIELDS = ("AroR", "HA", "MW", "HBA", "HBD", "TPSA")


def identity_key(row):
    label = str(row.get("compound_id") if row.get("compound_id") is not None else "").strip().casefold()
    if not label:
        return None
    values = []
    for key in STRUCTURE_FIELDS:
        value = row.get(key)
        if isinstance(value, bool) or value is None:
            return None
        try:
            value = float(value)
        except (TypeError, ValueError):
            return None
        if not math.isfinite(value):
            return None
        values.append(format(value, ".12g"))
    return hashlib.sha256(json.dumps([label, values], ensure_ascii=False).encode()).hexdigest()


@lru_cache(maxsize=2)
def _read_index(path, modified, size):
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("version") != 1:
        raise ValueError("Unsupported structure index version. Rebuild the index.")
    return data


def load_index():
    if not INDEX_PATH.exists():
        return {"version": 1, "records": {}}
    stat = INDEX_PATH.stat()
    return _read_index(str(INDEX_PATH), stat.st_mtime_ns, stat.st_size)


def lookup_structure(row, index=None):
    index = load_index() if index is None else index
    records = index.get("records", {})
    key = identity_key(row)
    if key:
        return records.get(key)
    # A name-only request is supported only when no descriptors were supplied
    # and every approved record for this name identifies the same molecular graph.
    if any(row.get(k) is not None and str(row.get(k)).strip() for k in STRUCTURE_FIELDS):
        return None
    label = str(row.get("compound_id") or "").strip().casefold()
    matches = [r for r in records.values() if r["compound_id"].strip().casefold() == label
               and r.get("allow_name_lookup", False)]
    if matches and len({r["canonical_smiles"] for r in matches}) == 1:
        return matches[0]
    return None


def index_status():
    records = load_index().get("records", {})
    return {"records": len(records), "names": sorted({r["compound_id"] for r in records.values()})}
