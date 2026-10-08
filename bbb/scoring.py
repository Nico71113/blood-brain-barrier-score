"""Gupta et al. (2019), author-supplied active Excel calculator.

DOI: 10.1021/acs.jmedchem.9b01220; supplement s001, cells C3:C11.
Coefficients, branch boundaries, and weights follow that calculator exactly.
The polynomial outputs are deliberately not clipped or rescaled.
"""
import math

METHOD = "Gupta 2019 · author Excel calculator"
SOURCE = "https://acs.figshare.com/articles/dataset/10046255"
FIELDS = ("AroR", "HA", "MW", "HBA", "HBD", "TPSA", "pKa")
COUNTS = {"AroR", "HA", "HBA", "HBD"}
WEIGHTS = {"AroR": 1, "HA": 1, "MWHBN": 1.5, "TPSA": 2, "pKa": 0.5}


def number(value, label):
    if value is None or isinstance(value, bool) or str(value).strip() == "":
        raise ValueError(f"{label} is required.")
    try:
        result = float(value)
    except (ValueError, TypeError):
        raise ValueError(f"{label} must be a single number.") from None
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite.")
    return result


def score_descriptors(values):
    d = {key: number(values.get(key), key) for key in FIELDS}
    for key in COUNTS:
        if d[key] < 0 or not d[key].is_integer():
            raise ValueError(f"{key} must be a non-negative integer.")
        d[key] = int(d[key])
    if d["HA"] < 1 or d["MW"] <= 0 or d["TPSA"] < 0:
        raise ValueError("HA and MW must be positive; TPSA must be non-negative.")
    ar, ha, mw, hba, hbd, tpsa, pka = (d[k] for k in FIELDS)
    try:
        x = (hba + hbd) / math.sqrt(mw)
        p = {
            "AroR": {0: .336376, 1: .816016, 2: 1., 3: .691115, 4: .199399}.get(ar, 0.),
            "HA": (.0000443*ha**3 - .004556*ha**2 + .12775*ha - .463)/.624231 if 5 < ha <= 45 else 0.,
            "MWHBN": (26.733*x**3 - 31.495*x**2 + 9.5202*x - .1358)/.72258 if .05 < x <= .45 else 0.,
            "TPSA": (-.0067*tpsa + .9598)/.9598 if 0 < tpsa <= 120 else 0.,
            "pKa": (.00045068*pka**4 - .016331*pka**3 + .18618*pka**2 - .71043*pka + .8579)/.597488 if 3 < pka <= 11 else 0.,
        }
    except OverflowError:
        raise ValueError("Input values are too large.") from None
    contributions = [
        {"name": key, "value": x if key == "MWHBN" else d[key],
         "p": val, "weight": WEIGHTS[key], "contribution": val*WEIGHTS[key]}
        for key, val in p.items()
    ]
    score = sum(c["contribution"] for c in contributions)
    if not math.isfinite(score):
        raise ValueError("Input is outside the supported range.")
    return {"bbb_score": score, "descriptors": d, "mwhbn": x,
            "contributions": contributions, "method": METHOD, "source": SOURCE}
