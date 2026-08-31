#!/usr/bin/env python3
"""
Score candidate G4-folding JSON files against a reference.

Usage:
    python score_folds.py reference.json cand1.json cand2.json ... \
        [--csv results.csv] [--filter-valid] [--sort {score,file}] [--desc]

    python find_candidates | python score_folds.py reference.json - \
        --csv results.csv

A candidate is VALID only if every quadruplex tetrads' {nt1..nt4} nucleotide
sets match the reference (as sets, order-insensitive, across all helices).
Valid folds are then scored 0-100 via weighted deep comparison of everything
under helices[].
"""

import argparse
import csv
import json
import math
import sys

# ---------------------------------------------------------------- weights

DEFAULT_WEIGHTS = {
    # per-quadruplex criteria
    "tetrads_match": 10,
    "gba_classification_match": 8,
    "onz_match": 8,
    "planarity_deviation": 6,
    "onzm": 4,
    "handedness": 6,
    "tetradPolarities": 5,
    "strandPolarities": 5,
    "loop_classification": 6,
    "quad_gba": 5,
    "tracts": 6,
    "path": 5,
    "bulges": 4,
    "loops": 6,
    "tetrad_pair_direction": 4,
    "tetrad_pair_geometry": 6,
    # document-level terms
    "dotBracket": 5,
    "quadruplexDotBracket": 5,
    # baseline ceiling (points) for an invalid-but-nonempty candidate
    "invalid_baseline": 15.0,
}


def load_weights(path):
    """Load weight overrides from a JSON file; missing keys use defaults."""
    with open(path) as f:
        overrides = json.load(f)
    unknown = set(overrides) - set(DEFAULT_WEIGHTS)
    if unknown:
        raise ValueError(f"Unknown weight keys: {sorted(unknown)}")
    w = dict(DEFAULT_WEIGHTS)
    w.update(overrides)
    return w


# ---------------------------------------------------------------- helpers


def tetrad_nts(tetrad):
    return frozenset(tetrad.get(k) for k in ("nt1", "nt2", "nt3", "nt4"))


def quad_tetrads(quad):
    return tuple(sorted(tetrad_nts(t) for t in quad.get("tetrads", [])))


def extract_quads(doc):
    """All quadruplexes across all helices."""
    quads = []
    helices = doc.get("helices", [])
    if isinstance(helices, dict):
        helices = [helices]
    for helix in helices:
        quads.extend(helix.get("quadruplexes", []))
    return quads


def sim(a, b):
    """Similarity of two numeric values relative to reference a (0..1)."""
    if a == 0:
        return 1.0 if b == 0 else 0.0
    return max(0.0, 1.0 - abs(a - b) / abs(a))


def rmsd(vals):
    if not vals:
        return 0.0
    return math.sqrt(sum(v * v for v in vals) / len(vals))


def _canon(x):
    """Canonical, comparable representation of dicts/lists/sets/scalars."""
    if isinstance(x, (set, frozenset)):
        return ("set", sorted(map(_canon, x)))
    if isinstance(x, dict):
        return ("dict", sorted((k, _canon(v)) for k, v in x.items()))
    if isinstance(x, (list, tuple)):
        return ("list", [_canon(v) for v in x])
    return x


def multiset_eq(a, b):
    return sorted(map(repr, map(_canon, a))) == sorted(map(repr, map(_canon, b)))


# ------------------------------------------------- matching quads by tetrads


def match_quads(ref_quads, cand_quads):
    """Greedy 1-1 matching by tetrad Jaccard overlap."""
    pairs = []
    used = set()
    for rq in ref_quads:
        best, best_j = None, -1.0
        for ci, cq in enumerate(cand_quads):
            if ci in used:
                continue
            rt = set(quad_tetrads(rq))
            ct = set(quad_tetrads(cq))
            j = len(rt & ct) / max(1, len(rt | ct))
            if j > best_j:
                best, best_j = ci, j
        pairs.append((rq, cand_quads[best] if best is not None else None))
        if best is not None:
            used.add(best)
    return pairs


# ------------------------------------------------------------ comparison


def compare_quad(rq, cq, r_tetrad_pairs=None, c_tetrad_pairs=None, weights=None):
    """Return (score_fraction 0..1, detail_dict) for one ref quad vs matched
    candidate quad. tetradPairs live on the helix, so they are passed in."""
    w = weights if weights is not None else DEFAULT_WEIGHTS
    detail = {}
    state = {"total": 0, "score": 0.0}
    flags = {}

    def add(name, weight, ok):
        state["total"] += weight
        state["score"] += weight * (1.0 if ok else 0.0)
        flags[name] = ok

    if cq is None:
        return 0.0, {"missing_quadruplex": True}

    # --- tetrads multiset (gate; also contributes to score)
    rt = [tetrad_nts(t) for t in rq.get("tetrads", [])]
    ct = [tetrad_nts(t) for t in cq.get("tetrads", [])]
    add("tetrads_match", w["tetrads_match"], multiset_eq(rt, ct))

    # --- per-tetrad fields
    r_by = {tetrad_nts(t): t for t in rq.get("tetrads", [])}
    c_by = {tetrad_nts(t): t for t in cq.get("tetrads", [])}
    gba_ok = onz_ok = True
    plan_deltas = []
    for t in rt:
        r, c = r_by.get(t), c_by.get(t)
        if r is None or c is None:
            gba_ok = onz_ok = False
            plan_deltas.append(1.0)
            continue
        gba_ok &= r.get("gbaClassification") == c.get("gbaClassification")
        onz_ok &= r.get("onz") == c.get("onz")
        plan_deltas.append(
            1.0
            - sim(r.get("planarityDeviation", 0.0), c.get("planarityDeviation", 0.0))
        )
    add("gba_classification_match", 8, gba_ok)
    add("onz_match", 8, onz_ok)
    pd = rmsd(plan_deltas)
    state["total"] += w["planarity_deviation"]
    state["score"] += w["planarity_deviation"] * (1.0 - pd)
    detail["planarity_delta"] = round(pd, 4)

    # --- quad-level scalars
    for name in ("onzm", "handedness", "tetradPolarities", "strandPolarities"):
        add(name + "_match", w[name], rq.get(name) == cq.get(name))

    # --- loop classification & gba classification
    lc = rq.get("loopClassification") or {}
    cc = cq.get("loopClassification") or {}
    add(
        "loop_classification_match",
        w["loop_classification"],
        lc.get("classification") == cc.get("classification")
        and lc.get("loopProgression") == cc.get("loopProgression"),
    )
    add(
        "quad_gba_match",
        w["quad_gba"],
        rq.get("gbaClassification") == cq.get("gbaClassification"),
    )

    # --- tracts: order within tract matters, list of tracts doesn't
    def tracts_key(q):
        return sorted(tuple(t) for t in q.get("tracts", []))

    add("tracts_match", w["tracts"], tracts_key(rq) == tracts_key(cq))

    # --- path (order-sensitive)
    add("path_match", w["path"], rq.get("path") == cq.get("path"))

    # --- bulges & loops as multisets
    add(
        "bulges_match",
        w["bulges"],
        multiset_eq(rq.get("bulges", []), cq.get("bulges", [])),
    )
    add(
        "loops_match", w["loops"], multiset_eq(rq.get("loops", []), cq.get("loops", []))
    )

    # --- tetradPairs (live on the helix; passed in from caller)
    rtp = r_tetrad_pairs
    ctp = c_tetrad_pairs
    if rtp is not None or ctp is not None:
        rtp = rtp or []
        ctp = ctp or []
        key = lambda p: (p.get("tetrad1"), p.get("tetrad2"))
        c_by_pair = {key(p): p for p in ctp}
        direction_ok = True
        rise_d, twist_d = [], []
        for p in rtp:
            c = c_by_pair.get(key(p)) or c_by_pair.get((key(p)[1], key(p)[0]))
            if c is None:
                direction_ok = False
                rise_d.append(1.0)
                twist_d.append(1.0)
                continue
            direction_ok &= p.get("direction") == c.get("direction")
            rise_d.append(1.0 - sim(p.get("rise", 0.0), c.get("rise", 0.0)))
            twist_d.append(1.0 - sim(p.get("twist", 0.0), c.get("twist", 0.0)))
        add("tetrad_pair_direction_match", w["tetrad_pair_direction"], direction_ok)
        gw = w["tetrad_pair_geometry"]
        state["total"] += gw
        state["score"] += gw * (
            (1.0 - rmsd(rise_d)) * 0.5 + (1.0 - rmsd(twist_d)) * 0.5
        )
        detail["rise_rmsd"] = round(rmsd(rise_d), 4)
        detail["twist_rmsd"] = round(rmsd(twist_d), 4)

    detail.update(flags)
    return (state["score"] / state["total"] if state["total"] else 1.0), detail


# ------------------------------------------------------------ main scoring


def score_candidate(ref, cand, weights=None):
    ref_quads = extract_quads(ref)
    cand_quads = extract_quads(cand)

    # ---- validity gate: same multiset of tetrad nt-sets across all quads
    ref_tetrads = sorted(t for q in ref_quads for t in quad_tetrads(q))
    cand_tetrads = sorted(t for q in cand_quads for t in quad_tetrads(q))
    valid = ref_tetrads == cand_tetrads

    w = weights if weights is not None else DEFAULT_WEIGHTS
    result = {
        "valid": valid,
        "score": None,
        "tetrads_match": valid,
        "details": {},
    }

    # zero quadruplexes in the candidate => hard zero
    if not cand_quads:
        result["score"] = 0.0
        result["details"]["no_quadruplexes"] = True
        return result

    # ---- tetradPairs per helix (matched by helix index)
    def helix_pairs(doc, idx):
        helices = doc.get("helices", [])
        if isinstance(helices, dict):
            helices = [helices]
        return helices[idx].get("tetradPairs") if idx < len(helices) else None

    # ---- deep comparison
    pairs = match_quads(ref_quads, cand_quads)
    wsum, ssum = 0, 0.0
    for i, (rq, cq) in enumerate(pairs):
        frac, detail = compare_quad(
            rq, cq, helix_pairs(ref, i), helix_pairs(cand, i), weights=w
        )
        weight = len(rq.get("tetrads", [])) or 1
        wsum += weight
        ssum += weight * frac
        for k, v in detail.items():
            result["details"][f"q{i}_{k}"] = v
    quad_score = 100.0 * ssum / wsum if wsum else 100.0

    # candidate has quads but none match the reference tetrads:
    # give a small baseline so it's not zero, scaled by how close
    # the quad count is to the reference
    if not valid:
        n_ref, n_cand = len(ref_quads), len(cand_quads)
        count_sim = 1.0 - abs(n_ref - n_cand) / max(n_ref, n_cand, 1)
        quad_score = min(quad_score, w["invalid_baseline"] * count_sim)

    # ---- dot brackets: configurable weights, dynamically normalized
    db_match = ref.get("dotBracket") == cand.get("dotBracket")
    qdb_match = ref.get("quadruplexDotBracket") == cand.get("quadruplexDotBracket")
    result["details"]["dot_bracket_match"] = db_match
    result["details"]["quad_dot_bracket_match"] = qdb_match

    # dynamic normalization: weighted sum over total weight => 0..100
    total_w = 100.0 + w["dotBracket"] + w["quadruplexDotBracket"]
    ssum = (
        100.0 * quad_score
        + w["dotBracket"] * 100.0 * db_match
        + w["quadruplexDotBracket"] * 100.0 * qdb_match
    )
    result["score"] = round(ssum / total_w, 2)
    return result


# ------------------------------------------------------------------ CLI

CSV_FIELDS = [
    "file",
    "valid",
    "score",
    "tetrads_match",
    "dot_bracket_match",
    "quad_dot_bracket_match",
]


def flatten(res):
    row = {
        "file": res["file"],
        "valid": str(res["valid"]),
        "score": res["score"],
        "tetrads_match": str(res["tetrads_match"]),
    }
    for k, v in res["details"].items():
        row[k] = str(v) if isinstance(v, bool) else v
    return row


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("reference", help="Reference JSON file (correct values)")
    ap.add_argument(
        "candidates",
        nargs="*",
        help="Candidate JSON files; use '-' to read paths from stdin",
    )
    ap.add_argument("--csv", help="Write comparison results to this CSV file")
    ap.add_argument(
        "--weights",
        help="JSON file with weight overrides (keys must match DEFAULT_WEIGHTS)",
    )
    ap.add_argument(
        "--filter-valid", action="store_true", help="Only include valid folds in output"
    )
    ap.add_argument("--sort", choices=["score", "file"], default="score")
    ap.add_argument(
        "--desc",
        action="store_true",
        help="Sort descending (score sorts descending by default)",
    )
    args = ap.parse_args()

    weights = load_weights(args.weights) if args.weights else None

    with open(args.reference) as f:
        ref = json.load(f)

    paths = list(args.candidates)
    if "-" in paths or (not paths and not sys.stdin.isatty()):
        paths = [p for p in paths if p != "-"] + [
            line.strip() for line in sys.stdin if line.strip()
        ]

    results = []
    for path in paths:
        try:
            with open(path) as f:
                cand = json.load(f)
            res = score_candidate(ref, cand, weights=weights)
        except Exception as e:
            res = {
                "file": path,
                "valid": False,
                "score": 0.0,
                "tetrads_match": False,
                "details": {"error": str(e)},
            }
        res["file"] = path
        results.append(res)

    if args.filter_valid:
        results = [r for r in results if r["valid"]]

    if args.sort == "score":
        results.sort(key=lambda r: (r["score"], r["file"]), reverse=not args.desc)
    else:
        results.sort(key=lambda r: r["file"], reverse=args.desc)

    # collect all detail keys so nothing is dropped in the CSV
    fields = list(CSV_FIELDS)
    for r in results:
        for k in r["details"]:
            if k not in fields:
                fields.append(k)

    rows = [flatten(r) for r in results]

    out = open(args.csv, "w", newline="") if args.csv else sys.stdout
    try:
        w = csv.DictWriter(out, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    finally:
        if args.csv:
            out.close()
            print(f"Wrote {len(rows)} rows to {args.csv}")


if __name__ == "__main__":
    main()
