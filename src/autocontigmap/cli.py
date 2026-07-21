#!/usr/bin/env python3
"""
autocontigmap — CLI for motif gap estimation.

Usage
-----
autocontigmap <pdb_file> <res_min> <res_max> [--pickle-file NAME_OR_PATH]

Arguments
---------
pdb_file      Path to the motif PDB file
res_min       Minimum total protein length in residues (integer, 10-499)
res_max       Maximum total protein length in residues (integer, 10-499)

Output (stdout)
---------------
Standard (default) output, comma-separated, with fixed motif/target chain
spans included alongside the estimated gap-fill ranges:

    "contig": "8-15,A18-25,16-30,A47-54,16-30,A92-99"
    "length": "150-200"

Pass --rfd1 to instead print the equivalent old-style RFdiffusion(1)
contig, slash-separated:

    contigmap.contigs=[8-15/A18-25/16-30/A47-54/16-30/A92-99]
    contigmap.length="150-200"

"length"/contigmap.length is just res_min-res_max echoed back, for pinning
the overall design length directly (the default terminal budget below
otherwise only pins a 0-upper-bound range, not an exact total).

Terminal (N-/C-terminal tail) budget defaults to the simple 0-term_hi form
(term_hi = floor((res_max - motif - gap_lo_sum) / 2)). Pass
--strict-terminals for the old res_min/res_max conflict-clamped form
instead (clamped to a single value at both ends if the two conflict, with a
warning printed to stderr). A terminal segment that computes to exactly 0-0
is omitted from the contig entirely rather than printed as a literal "0".

If a gap's estimated Cα-Cα distance is beyond the 95th percentile of gap
sizes seen in the checkpoint data for this res_min-res_max range, a warning
is printed to stderr -- the estimate is based on thin data out there, and a
larger design is probably necessary. If the requested res_min/res_max range
can't fit the motif and its gaps at all, that's a hard error (or, if it's
merely too tight, a warning) -- see "Length-range validation" below.

Diagnostics (stderr)
---------------------
Every run starts with what it found ("Chains found: ...", then how the
motif was identified), before any gap processing. All diagnostic lines use
the same two-space indent and "Warning: ..." phrasing; a blank line follows
the per-gap listing before any length-range warning. Hard errors are
checked first, before any gap-related warning would print.

Length-range validation
------------------------
Hard error if res_min is smaller than the motif's own residue count
(excluding gaps) -- checked before any gap processing, so nothing else
prints beforehand. Warning (not fatal) if res_min clears that bar but
res_max is still too tight for even the smallest per-gap estimate -- the
terminal budget clamps to 0 and the contig is still printed.

Chain selection
---------------
Any chain that contains at least one internal gap (chain break) is treated
as a motif chain and gap-filled. If no chain has a gap and --chain-order
isn't given, this is an error -- there is nothing to scaffold.

Default mode (no --chain-order): every chain with an internal gap is
gap-filled independently (its own terminal budget, its own designed
segment); chains without a gap are carried through unchanged as fixed
spans. All segments are joined with a hard chain break (/0 ) in the
output's chain order.

If the motif's segments are instead split across separate PDB chains (one
chain per segment, no internal gaps), pass --chain-order with the segment
chains in contig order (e.g. --chain-order B,A). Gaps are then measured
between the last residue of each chain and the first residue of the next,
and those chains are scaffolded into a single continuous designed chain
(no break between them) sharing one terminal budget. Any other chains
present are still carried through unchanged as fixed spans, each behind its
own chain break.
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import Selection

from .core import DEFAULT_CHECKPOINT, aggregate_by_residue_range, ca_distances, load_pdb, load_pickle

GAP_SIZE_WARN_PERCENTILE = 0.95


class AutoContigmapError(Exception):
    """Base for user-facing errors: caught in main(), printed to stderr, exit 1."""


class MotifNotFoundError(AutoContigmapError):
    pass


class LengthRangeError(AutoContigmapError):
    pass


def find_gaps_in_chain(chain):
    """Return list of (i_prev, i_curr, prev_resnum, curr_resnum) for each gap."""
    res_list = Selection.unfold_entities(chain, "R")
    gaps = []
    for i in range(1, len(res_list)):
        prev_num = res_list[i - 1].get_id()[1]
        curr_num = res_list[i].get_id()[1]
        if curr_num - prev_num > 1:
            gaps.append((i - 1, i, prev_num, curr_num))
    return gaps


def find_gapped_chains(structure):
    """[(chain, gaps), ...] for every chain in the structure that has at least one internal gap."""
    result = []
    for chain in structure[0].get_list():
        gaps = find_gaps_in_chain(chain)
        if gaps:
            result.append((chain, gaps))
    return result


def get_chain(structure, chain_id):
    try:
        return structure[0][chain_id]
    except KeyError:
        available = [c.id for c in structure[0]]
        raise ValueError(f"Chain '{chain_id}' not found in structure (available: {available})")


def chain_span(chain):
    """(first_resnum, last_resnum) of a chain."""
    res_list = Selection.unfold_entities(chain, "R")
    return res_list[0].get_id()[1], res_list[-1].get_id()[1]


def interchain_gap_distance(chain_prev, chain_curr):
    """Ca-Ca distance (Angstrom) between the last residue of chain_prev and the first residue of chain_curr."""
    prev_res = Selection.unfold_entities(chain_prev, "R")[-1]
    curr_res = Selection.unfold_entities(chain_curr, "R")[0]
    return np.linalg.norm(prev_res["CA"].coord - curr_res["CA"].coord)


def term_token(lo, hi):
    """Terminal segment string, or None if it's an empty (0-0) budget -- omit those entirely."""
    if lo == 0 and hi == 0:
        return None
    return str(lo) if lo == hi else f"{lo}-{hi}"


def compute_terminals(motif_residues, gap_lo_sum, gap_hi_sum, res_min, res_max):
    """
    Returns (term_lo, term_hi, clamped: bool).

    Normal case:
        term_hi = floor((res_max - motif - gap_lo_sum) / 2)
        term_lo = floor((res_min - motif - gap_hi_sum) / 2)

    Conflict: term_lo > term_hi or either is negative.
        Clamp both to floor((res_max - motif - gap_hi_sum) / 2).
        If even that is negative, clamp to 0.
    """
    term_hi = math.floor((res_max - motif_residues - gap_lo_sum) / 2)
    term_lo = math.floor((res_min - motif_residues - gap_hi_sum) / 2)

    if term_lo < 0 or term_lo > term_hi:
        fallback = math.floor((res_max - motif_residues - gap_hi_sum) / 2)
        fallback = max(0, fallback)
        return fallback, fallback, True

    return term_lo, term_hi, False


def compute_terminals_simple(motif_residues, gap_lo_sum, res_max):
    """
    Simplified terminal budget (the default): always 0-term_hi at both ends,
    skipping the res_min/res_max conflict clamping in compute_terminals().
    Meant to be paired with the "length" field pinning the exact overall
    scaffold length range directly, instead of relying on the terminal
    range to do it.

    term_hi = floor((res_max - motif - gap_lo_sum) / 2), floored at 0.
    """
    term_hi = math.floor((res_max - motif_residues - gap_lo_sum) / 2)
    return 0, max(0, term_hi)


def check_motif_fits(motif_residues, res_min):
    """Hard error -- checked before any gap processing/warnings, not just alongside them."""
    if res_min < motif_residues:
        raise LengthRangeError(
            f"res_min ({res_min}) is smaller than the motif's own residue count "
            f"({motif_residues}) -- the requested length range can't even fit the fixed "
            f"motif residues, let alone the gaps."
        )


def terminal_budget(motif_residues, gap_lo_sum, gap_hi_sum, res_min, res_max, use_strict):
    """Shared term_lo/term_hi resolution, printing a warning if the requested range is too tight."""
    if use_strict:
        term_lo, term_hi, clamped = compute_terminals(motif_residues, gap_lo_sum, gap_hi_sum, res_min, res_max)
        if clamped:
            suggested_max = motif_residues + gap_hi_sum + 2 * max(term_lo, 1)
            print(
                f"  Warning: gap estimates ({gap_lo_sum}-{gap_hi_sum} aa) leave no consistent "
                f"terminal budget within {res_min}-{res_max} aa; terminals clamped to {term_lo} "
                f"each -- try res_max >= {suggested_max}.",
                file=sys.stderr,
            )
    else:
        term_lo, term_hi = compute_terminals_simple(motif_residues, gap_lo_sum, res_max)
        if res_max - motif_residues - gap_lo_sum < 0:
            suggested_max = motif_residues + gap_lo_sum
            print(
                f"  Warning: even the smallest gap estimates ({gap_lo_sum} aa total) push the "
                f"minimum feasible length to {suggested_max} aa, above res_max ({res_max}) -- "
                f"a bigger design length range is probably necessary (res_max >= {suggested_max}).",
                file=sys.stderr,
            )

    return term_lo, term_hi


def gap_size_warn_threshold(agg):
    """Gap size (Å) at the GAP_SIZE_WARN_PERCENTILE of the aggregated gap-size distribution, or None."""
    a = np.cumsum(agg["counts"])
    norm = (a - a.min()) / (a.max() - a.min())
    idx = np.where(norm > GAP_SIZE_WARN_PERCENTILE)[0]
    return int(idx[0]) + 1 if len(idx) else None


def warn_if_gap_too_large(label, gap_ang, threshold, res_min, res_max):
    if threshold is not None and gap_ang > threshold:
        pct = int(GAP_SIZE_WARN_PERCENTILE * 100)
        print(
            f"  Warning: gap {label} is {gap_ang} Å, beyond the {pct}th percentile "
            f"({threshold} Å) of gap sizes seen in our data for {res_min}-{res_max}-residue "
            f"proteins -- a bigger design length range is probably necessary.",
            file=sys.stderr,
        )


# ---------------------------------------------------------------------------
# Contig assembly, shared between the standard and --rfd1 output styles
# (they differ only in the intra-chain separator: "," vs "/")
# ---------------------------------------------------------------------------

def _gapped_chain_segment(chain, gaps, res_min, res_max, gap_size_data, use_strict_terminals, sep):
    """One chain with internal gaps: [term/]fixed/gap/fixed[/.../term], joined by sep."""
    cid = chain.get_id()
    first_res, last_res = chain_span(chain)
    motif_residues = len(Selection.unfold_entities(chain, "R"))
    check_motif_fits(motif_residues, res_min)

    distances = ca_distances(chain)
    agg = aggregate_by_residue_range(gap_size_data, res_min, res_max)
    threshold = gap_size_warn_threshold(agg)

    gap_estimates = []  # (aa_low, aa_high, prev_resnum, curr_resnum)
    for (i_prev, i_curr, prev_resnum, curr_resnum) in gaps:
        gap_ang = math.floor(distances[i_prev][i_curr])
        aa_low = int(np.floor(agg["Q0.4"][gap_ang - 1]))
        aa_high = int(np.floor(agg["Q0.6"][gap_ang - 1]))
        gap_estimates.append((aa_low, aa_high, prev_resnum, curr_resnum))
        print(f"  Gap {cid}{prev_resnum}-{curr_resnum}: {gap_ang} Å  ->  {aa_low}-{aa_high} residues", file=sys.stderr)
        warn_if_gap_too_large(f"{cid}{prev_resnum}-{curr_resnum}", gap_ang, threshold, res_min, res_max)
    print(file=sys.stderr)

    gap_lo_sum = sum(lo for lo, hi, *_ in gap_estimates)
    gap_hi_sum = sum(hi for lo, hi, *_ in gap_estimates)
    term_lo, term_hi = terminal_budget(
        motif_residues, gap_lo_sum, gap_hi_sum, res_min, res_max, use_strict_terminals
    )

    parts = []
    lead = term_token(term_lo, term_hi)
    if lead is not None:
        parts.append(lead)
    cursor = first_res
    for lo, hi, prev_resnum, curr_resnum in gap_estimates:
        parts.append(f"{cid}{cursor}-{prev_resnum}")
        parts.append(f"{lo}-{hi}")
        cursor = curr_resnum
    parts.append(f"{cid}{cursor}-{last_res}")
    trail = term_token(term_lo, term_hi)
    if trail is not None:
        parts.append(trail)
    return sep.join(parts)


def _chain_order_segment(chains, chain_ids, res_min, res_max, gap_size_data, use_strict_terminals, sep):
    """Chains merged in chain_ids order via inter-chain gaps: [term/]fixed/gap/fixed[/.../term], joined by sep."""
    motif_residues = sum(len(Selection.unfold_entities(c, "R")) for c in chains)
    check_motif_fits(motif_residues, res_min)

    agg = aggregate_by_residue_range(gap_size_data, res_min, res_max)
    threshold = gap_size_warn_threshold(agg)

    gap_estimates = []  # (aa_low, aa_high)
    for i in range(len(chains) - 1):
        gap_ang = math.floor(interchain_gap_distance(chains[i], chains[i + 1]))
        aa_low = int(np.floor(agg["Q0.4"][gap_ang - 1]))
        aa_high = int(np.floor(agg["Q0.6"][gap_ang - 1]))
        gap_estimates.append((aa_low, aa_high))
        label = f"{chain_ids[i]}(last)-{chain_ids[i + 1]}(first)"
        print(f"  Gap {label}: {gap_ang} Å  ->  {aa_low}-{aa_high} residues", file=sys.stderr)
        warn_if_gap_too_large(label, gap_ang, threshold, res_min, res_max)
    print(file=sys.stderr)

    gap_lo_sum = sum(lo for lo, hi in gap_estimates)
    gap_hi_sum = sum(hi for lo, hi in gap_estimates)
    term_lo, term_hi = terminal_budget(
        motif_residues, gap_lo_sum, gap_hi_sum, res_min, res_max, use_strict_terminals
    )

    parts = []
    lead = term_token(term_lo, term_hi)
    if lead is not None:
        parts.append(lead)
    for i, chain in enumerate(chains):
        first_res, last_res = chain_span(chain)
        parts.append(f"{chain.get_id()}{first_res}-{last_res}")
        if i < len(chains) - 1:
            parts.append(f"{gap_estimates[i][0]}-{gap_estimates[i][1]}")
    trail = term_token(term_lo, term_hi)
    if trail is not None:
        parts.append(trail)
    return sep.join(parts)


def build_contig(structure, res_min, res_max, gap_size_data, chain_order=None, use_strict_terminals=False, sep=","):
    """
    Assemble the full contig body (without the surrounding "contigmap.contigs=[...]"
    or '"contig": "..."' wrapper): one designed/fixed segment per output
    chain, joined by a hard chain break ("/0 "). `sep` is the intra-chain
    separator: "," for the standard style, "/" for --rfd1.
    """
    all_chains = list(structure[0].get_list())
    print(f"  Chains found: {', '.join(c.get_id() for c in all_chains)}", file=sys.stderr)

    if chain_order:
        chain_ids = [c.strip() for c in chain_order.split(",") if c.strip()]
        if not chain_ids:
            raise ValueError("--chain-order must list at least one chain.")
        print(
            f"  Motif identified via --chain-order: {', '.join(chain_ids)} "
            f"(merged into one chain, gaps measured between consecutive chains)",
            file=sys.stderr,
        )
        ordered_chains = [get_chain(structure, cid) for cid in chain_ids]

        motif_segment = _chain_order_segment(
            ordered_chains, chain_ids, res_min, res_max, gap_size_data, use_strict_terminals, sep
        )

        used_ids = set(chain_ids)
        groups, inserted = [], False
        for chain in all_chains:
            if chain.get_id() in used_ids:
                if not inserted:
                    groups.append(motif_segment)
                    inserted = True
                continue
            first_res, last_res = chain_span(chain)
            groups.append(f"{chain.get_id()}{first_res}-{last_res}")
        return "/0 ".join(groups)

    gapped = find_gapped_chains(structure)
    if not gapped:
        raise MotifNotFoundError(
            "Not motif chain identified, motif chains have to have gaps. "
            "If the motif is defined over multiple chains use --chain-order"
        )
    gaps_by_id = {chain.get_id(): gaps for chain, gaps in gapped}
    print(
        f"  Motif identified automatically (has internal gap(s)): {', '.join(gaps_by_id)}",
        file=sys.stderr,
    )

    groups = []
    for chain in all_chains:
        cid = chain.get_id()
        if cid in gaps_by_id:
            groups.append(
                _gapped_chain_segment(
                    chain, gaps_by_id[cid], res_min, res_max, gap_size_data, use_strict_terminals, sep
                )
            )
        else:
            first_res, last_res = chain_span(chain)
            groups.append(f"{cid}{first_res}-{last_res}")
    return "/0 ".join(groups)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Estimate gap sizes for a motif PDB and output an RFdiffusion "
            "contig string including N/C-terminal tail budgets."
        )
    )
    parser.add_argument("pdb_file", help="Path to the motif PDB file")
    parser.add_argument("res_min", type=int, help="Minimum total protein length (10-499)")
    parser.add_argument("res_max", type=int, help="Maximum total protein length (10-499)")
    parser.add_argument(
        "--pickle-file",
        default=DEFAULT_CHECKPOINT,
        help=(
            "Checkpoint dataset to use: one of the bundled variants "
            "(gyr [default], gyr_ss_2, standard, surface) or a path to an "
            "external .pkl file."
        ),
    )
    parser.add_argument(
        "--rfd1", action="store_true",
        help=(
            "Print the old-style RFdiffusion(1) contig "
            "('contigmap.contigs=[...]'/'contigmap.length=...') instead of "
            "the standard '\"contig\": ...'/'\"length\": ...' output."
        ),
    )
    parser.add_argument(
        "--strict-terminals", action="store_true",
        help=(
            "Use the res_min/res_max conflict-clamped terminal budget instead of the "
            "default simple 0-term_hi one (term_hi from the res_max/gap_lo_sum formula)."
        ),
    )
    parser.add_argument(
        "--chain-order",
        help=(
            "Comma-separated chain IDs in contig order (e.g. 'B,A'), for motifs whose "
            "segments are split across separate PDB chains. Gaps are measured between "
            "the last residue of each chain and the first residue of the next."
        ),
    )
    args = parser.parse_args()

    if not Path(args.pdb_file).exists():
        print(f"Error: pdb_file not found: {args.pdb_file}", file=sys.stderr)
        sys.exit(1)

    if not (10 <= args.res_min <= 499 and 10 <= args.res_max <= 499):
        print("Error: res_min and res_max must both be in [10, 499].", file=sys.stderr)
        sys.exit(1)

    if args.res_min > args.res_max:
        print("Error: res_min must be <= res_max.", file=sys.stderr)
        sys.exit(1)

    gap_size_data = load_pickle(args.pickle_file)
    structure = load_pdb(args.pdb_file)

    try:
        body = build_contig(
            structure, args.res_min, args.res_max, gap_size_data,
            chain_order=args.chain_order, use_strict_terminals=args.strict_terminals,
            sep="/" if args.rfd1 else ",",
        )
    except AutoContigmapError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    length_range = f"{args.res_min}-{args.res_max}"
    if args.rfd1:
        print(f"contigmap.contigs=[{body}]")
        print(f'contigmap.length="{length_range}"')
    else:
        print(f'"contig": "{body}"')
        print(f'"length": "{length_range}"')


if __name__ == "__main__":
    main()
