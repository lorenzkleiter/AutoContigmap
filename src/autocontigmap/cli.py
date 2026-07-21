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
    contig="5-15,60-100,5-15"

where the first and last entries are the N- and C-terminal tail budgets
and the middle entries are the estimated residue counts per internal gap.

If the gap estimates leave no room for terminals within res_max, the
terminals are clamped to floor((res_max - motif - sum(gap_hi)) / 2) and
a warning is printed to stderr.

Chain selection
---------------
By default, the motif chain (the one to be scaffolded) is chosen automatically:
  1. Prefer the chain that contains at least one internal gap (chain break).
  2. If no chain has a gap, or multiple chains tie, pick the shortest one.

If the motif's segments are instead split across separate PDB chains (one
chain per segment), pass --chain-order with the segment chains in contig
order (e.g. --chain-order B,A). Gaps are then measured between the last
residue of each chain and the first residue of the next, instead of within
a single chain.
"""

import argparse
import math
import sys
from pathlib import Path

import numpy as np
from Bio.PDB import Selection

from .core import DEFAULT_CHECKPOINT, aggregate_by_residue_range, ca_distances, load_pdb, load_pickle


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


def get_chain(structure, chain_id):
    try:
        return structure[0][chain_id]
    except KeyError:
        available = [c.id for c in structure[0]]
        raise ValueError(f"Chain '{chain_id}' not found in structure (available: {available})")


def interchain_gap_distance(chain_prev, chain_curr):
    """Ca-Ca distance (Angstrom) between the last residue of chain_prev and the first residue of chain_curr."""
    prev_res = Selection.unfold_entities(chain_prev, "R")[-1]
    curr_res = Selection.unfold_entities(chain_curr, "R")[0]
    return np.linalg.norm(prev_res["CA"].coord - curr_res["CA"].coord)


def select_motif_chain(structure):
    """
    Pick the motif chain to scaffold:
      1. Prefer chains with at least one internal gap.
      2. Among tied candidates (all have gaps, or none do), pick the shortest.
    Returns (chain, gaps) where gaps is the list from find_gaps_in_chain().
    """
    chains = list(structure[0].get_list())

    chains_with_gaps = []
    chains_without = []
    for chain in chains:
        gaps = find_gaps_in_chain(chain)
        if gaps:
            chains_with_gaps.append((chain, gaps))
        else:
            chains_without.append((chain, gaps))

    candidates = chains_with_gaps if chains_with_gaps else chains_without
    candidates.sort(key=lambda cg: len(Selection.unfold_entities(cg[0], "R")))
    chosen_chain, chosen_gaps = candidates[0]

    if len(chains) > 1:
        print(
            f"  Motif chain selected: {chosen_chain.get_id()} "
            f"({'has gap(s)' if chosen_gaps else 'no gaps, shortest chain'}, "
            f"{len(Selection.unfold_entities(chosen_chain, 'R'))} residues)",
            file=sys.stderr,
        )

    return chosen_chain, chosen_gaps


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
    Simplified terminal budget: always 0-term_hi at both ends, skipping the
    res_min/res_max conflict clamping in compute_terminals(). Intended for
    callers that pin the exact scaffold length via contigmap.length instead
    of relying on the terminal range, so res_min == res_max anyway.

    term_hi = floor((res_max - motif - gap_lo_sum) / 2), floored at 0.
    """
    term_hi = math.floor((res_max - motif_residues - gap_lo_sum) / 2)
    return 0, max(0, term_hi)


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
        "--simple-terminals", action="store_true",
        help=(
            "Skip the res_min/res_max conflict clamping and always emit 0-term_hi "
            "at both ends (term_hi from the res_max/gap_lo_sum formula)."
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

    gap_estimates = []  # list of (aa_low, aa_high, gap_ang, prev_label, curr_label)

    if args.chain_order:
        chain_ids = [c.strip() for c in args.chain_order.split(",") if c.strip()]
        if len(chain_ids) < 1:
            print("Error: --chain-order must list at least one chain.", file=sys.stderr)
            sys.exit(1)
        chains = [get_chain(structure, cid) for cid in chain_ids]
        motif_residues = sum(len(Selection.unfold_entities(c, "R")) for c in chains)

        print(f"  Motif chains (contig order): {','.join(chain_ids)}", file=sys.stderr)
        print(f"  Motif residues (excluding gaps): {motif_residues}", file=sys.stderr)

        agg = aggregate_by_residue_range(gap_size_data, args.res_min, args.res_max)
        for i in range(len(chains) - 1):
            gap_ang = math.floor(interchain_gap_distance(chains[i], chains[i + 1]))
            aa_low = int(np.floor(agg["Q0.4"][gap_ang - 1]))
            aa_high = int(np.floor(agg["Q0.6"][gap_ang - 1]))
            gap_estimates.append((aa_low, aa_high, gap_ang, chain_ids[i], chain_ids[i + 1]))
            print(
                f"  Gap {chain_ids[i]}(last)-{chain_ids[i + 1]}(first): {gap_ang} Å  ->  "
                f"{aa_low}-{aa_high} residues",
                file=sys.stderr,
            )
    else:
        motif_chain, raw_gaps = select_motif_chain(structure)
        res_list = Selection.unfold_entities(motif_chain, "R")
        motif_residues = len(res_list)

        print(f"  Motif residues (excluding gaps): {motif_residues}", file=sys.stderr)

        for (i_prev, i_curr, prev_resnum, curr_resnum) in raw_gaps:
            distances = ca_distances(motif_chain)
            gap_ang = math.floor(distances[i_prev][i_curr])
            agg = aggregate_by_residue_range(gap_size_data, args.res_min, args.res_max)
            aa_low = int(np.floor(agg["Q0.4"][gap_ang - 1]))
            aa_high = int(np.floor(agg["Q0.6"][gap_ang - 1]))
            gap_estimates.append((aa_low, aa_high, gap_ang, prev_resnum, curr_resnum))
            print(
                f"  Gap {prev_resnum}-{curr_resnum}: {gap_ang} Å  ->  {aa_low}-{aa_high} residues",
                file=sys.stderr,
            )

    gap_lo_sum = sum(lo for lo, hi, *_ in gap_estimates)
    gap_hi_sum = sum(hi for lo, hi, *_ in gap_estimates)

    if args.simple_terminals:
        term_lo, term_hi = compute_terminals_simple(motif_residues, gap_lo_sum, args.res_max)
        clamped = False
    else:
        term_lo, term_hi, clamped = compute_terminals(
            motif_residues, gap_lo_sum, gap_hi_sum, args.res_min, args.res_max
        )

    if clamped:
        suggested_max = motif_residues + gap_hi_sum + 2 * max(term_lo, 1)
        print(
            f"\n  WARNING: gap estimates ({gap_lo_sum}-{gap_hi_sum} aa) leave no consistent "
            f"terminal budget within the requested length range "
            f"({args.res_min}-{args.res_max} aa).\n"
            f"  Terminals clamped to {term_lo} residues each.\n"
            f"  To fully satisfy the length range, set res_max >= {suggested_max}.\n"
            f"  Gap ranges are unchanged.",
            file=sys.stderr,
        )

    def term_str(lo, hi):
        return str(lo) if lo == hi else f"{lo}-{hi}"

    parts = [term_str(term_lo, term_hi)]
    for lo, hi, *_ in gap_estimates:
        parts.append(f"{lo}-{hi}")
    parts.append(term_str(term_lo, term_hi))

    print(f'contig="{",".join(parts)}"')


if __name__ == "__main__":
    main()
