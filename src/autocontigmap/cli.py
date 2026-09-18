#!/usr/bin/env python3
"""
autocontigmap — CLI for motif gap estimation.

Usage
-----
autocontigmap <pdb_file> <res_min> [res_max] [--pickle-file NAME_OR_PATH]

Arguments
---------
pdb_file      Path to the motif PDB file
res_min       Minimum total protein length in residues (integer, 10-499)
res_max       Maximum total protein length in residues (integer, 10-499).
              If omitted, res_min is treated as a single fixed target length
              (internally res_max = res_min) and "length" is printed as that
              single value ("100") rather than a "min-max" range.

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

Terminal (N-/C-terminal tail) budget comes from the pooled permissive budget

    aa_term_low  = res_min - motif - gap_hi_sum
    aa_term_high = res_max - motif - gap_lo_sum

split over 2 * (number of designed segments) terminals, both bounds rounded
up. The lower bound is kept rather than discarded: it is pinned to 0 only
when aa_term_low < 0 while aa_term_low + sigma_gap >= 0, i.e. when every
requested length is still reachable with the gaps near their lower
estimates. If aa_term_low + sigma_gap < 0, res_min itself is unreachable and
that is a hard error. A terminal segment that computes to exactly 0-0 is omitted from the contig entirely
rather than printed as a literal "0".

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
prints beforehand. After gap estimation, two further hard errors: res_max
below motif + gap_lo_sum (no admissible configuration exists at all), and
res_min below motif + gap_lo_sum (the lower end of the requested range
cannot be reached at any terminal length). A negative lower terminal budget
that still leaves every length reachable is a warning, not an error.

Chain selection
---------------
Any chain that contains at least one internal gap (chain break) is treated
as a motif chain and gap-filled. If no chain has a gap and --chain-order
isn't given, this is an error -- there is nothing to scaffold.

Default mode (no --chain-order): every chain with an internal gap becomes
its own designed segment; chains without a gap are carried through
unchanged as fixed spans. All segments are joined with a hard chain break
(/0 ) in the output's chain order.

The terminal budget is NOT per segment. res_min/res_max are whole-design
lengths, so the motif count and the gap sums are pooled over every designed
segment and the resulting budget is split over all 2 * n_segments terminals.
Budgeting per segment would credit each one with the residues of all the
others: with two gapped 40 aa chains and one 85-125 aa gap each, requesting
300-400 used to emit 34-69 terminals, which can only express design totals
of 386-606 -- unable to reach 300 at all, and overshooting 400 by 200.

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

from .core import (
    DEFAULT_CHECKPOINT,
    MAX_CHECKPOINT_RESIDUES,
    GAP_SIZE_ERROR_PERCENTILE,
    GAP_SIZE_WARN_PERCENTILE,
    GapEstimateUnavailable,
    aggregate_by_residue_range,
    ca_distances,
    find_gaps_in_chain,
    check_gap_against_lengths,
    length_distance_thresholds,
    load_pdb,
    load_pickle,
    lookup_gap_estimate,
    pooled_terminal_budget,
    smallest_feasible_res_min,
    split_terminal_budget,
)


class AutoContigmapError(Exception):
    """Base for user-facing errors: caught in main(), printed to stderr, exit 1."""


class MotifNotFoundError(AutoContigmapError):
    pass


class LengthRangeError(AutoContigmapError):
    pass


class GapDataError(AutoContigmapError):
    pass


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


def _range_advice(minimal, recommended):
    """
    Closing clause of a length-range error: the smallest res_min that works, and
    the larger one that also clears the high end of the gap estimates.

    Both come from a re-estimating scan (see core.smallest_feasible_res_min),
    not from the floor computed for the window that just failed -- that floor
    moves as soon as res_min does, so quoting it back is bad advice.
    """
    if minimal is None:
        return " -- res_min and res_max both have to be increased."
    if recommended is None or recommended <= minimal:
        return f" -- try res_min >= {minimal}."
    return f" -- try res_min >= {minimal}, or >= {recommended} (recommended)."


def resolve_terminals(motif_residues, gap_lo_sum, gap_hi_sum, res_min, res_max,
                      n_segments=1, gap_size_data=None, gap_angstroms=()):
    """
    Per-terminal (term_lo, term_hi) from the pooled permissive budget in core.

    The lower bound is kept, not discarded: it is pinned to 0 only when
    aa_term_low < 0 while aa_term_low + sigma_gap >= 0, i.e. when every length
    in [res_min, res_max] is still reachable provided the gaps sample near
    their lower estimates. If aa_term_low + sigma_gap < 0, res_min itself is
    unreachable at any terminal assignment and that is a hard error rather than
    a clamp.

    Note that the emitted ranges still do not enforce [res_min, res_max] on
    their own -- the separate "length" field does that.
    """
    term_low, term_high, status = pooled_terminal_budget(
        motif_residues, gap_lo_sum, gap_hi_sum, res_min, res_max
    )
    if status in ("infeasible", "unreachable_min"):
        # One problem, not two. res_min <= res_max, so a res_max below the floor
        # ("infeasible") puts res_min below it too, and clearing res_min makes
        # res_max >= res_min >= floor automatic. Only res_min is worth scanning:
        # advising a bigger res_max would just hand back a range that fails the
        # other check on the next run.
        advice = _range_advice(
            smallest_feasible_res_min(
                gap_size_data, motif_residues, gap_angstroms, res_max
            ),
            smallest_feasible_res_min(
                gap_size_data, motif_residues, gap_angstroms, res_max, use_high=True
            ),
        )
        raise LengthRangeError(
            f"res_min ({res_min}) is below the shortest buildable design "
            f"({motif_residues} aa motif + {gap_lo_sum}-{gap_hi_sum} aa of gaps)"
            f"{advice}"
        )
    if status == "clamped":
        # The res_min that would make this "ok" is the same quantity the error
        # messages call "recommended": the smallest one whose own gap estimates
        # still fit at their high end. None when no res_min up to res_max does.
        comfortable = (
            None if gap_size_data is None
            else smallest_feasible_res_min(
                gap_size_data, motif_residues, gap_angstroms, res_max, use_high=True
            )
        )
        tip = "" if comfortable is None else f" Use res_min >= {comfortable} to avoid this."
        print(
            f"  Warning: the lower terminal budget is negative and was pinned to 0; "
            f"res_min ({res_min}) is still reachable, but only when the gaps sample near "
            f"their lower estimates ({gap_lo_sum} aa total).{tip}",
            file=sys.stderr,
        )

    return split_terminal_budget(term_low, term_high, n_segments)


def check_motif_fits(motif_residues, res_min):
    """Hard error -- checked before any gap processing/warnings, not just alongside them."""
    if res_min < motif_residues:
        raise LengthRangeError(
            f"res_min ({res_min}) is smaller than the motif's own residue count "
            f"({motif_residues}) -- the requested length range can't even fit the fixed "
            f"motif residues, let alone the gaps."
        )


def gap_size_thresholds(gap_size_data):
    """
    Per-length Cα-Cα distance thresholds at the warn and error percentiles.

    Computed from the raw checkpoint, before any aggregation over the requested
    design-length window, so the short end of the range is judged on its own
    data rather than on a pooled distribution that the longest chains dominate.
    """
    return (
        length_distance_thresholds(gap_size_data, GAP_SIZE_WARN_PERCENTILE),
        length_distance_thresholds(gap_size_data, GAP_SIZE_ERROR_PERCENTILE),
    )


def _verdict_message(label, gap_ang, verdict, res_min, res_max):
    pct = int(verdict.percentile * 100)
    lo, hi = verdict.bad_lengths
    span = f"{lo}" if lo == hi else f"{lo}-{hi}"
    msg = (
        f"gap {label} is {gap_ang} \N{ANGSTROM SIGN}, beyond the {pct}th percentile of "
        f"Cα-Cα distances for designs of {span} residues "
        f"({verdict.n_bad} of the {verdict.n_total} lengths in {res_min}-{res_max})"
    )
    if verdict.min_viable is None:
        return msg + (
            f". No chain length up to {MAX_CHECKPOINT_RESIDUES} supports a gap this wide; "
            f"the motif segments may be further apart than a single domain can span."
        )
    return msg + f". The shortest design length that supports it is {verdict.min_viable}."


def check_gap_size(label, gap_ang, thresholds, res_min, res_max):
    """
    Check one gap against every length in the window, before aggregation.

    Raises GapDataError past GAP_SIZE_ERROR_PERCENTILE, warns past
    GAP_SIZE_WARN_PERCENTILE, and names the lengths responsible plus the res_min
    that would clear the check.
    """
    verdict = check_gap_against_lengths(gap_ang, res_min, res_max, *thresholds)
    if verdict is None:
        return
    message = _verdict_message(label, gap_ang, verdict, res_min, res_max)
    if verdict.level == "error":
        raise GapDataError(message)
    print(f"  Warning: {message}", file=sys.stderr)


def _estimate(agg, gap_ang, label, res_min, res_max, thresholds):
    """Check the distance is plausible for every length in the window, then look it up."""
    check_gap_size(label, gap_ang, thresholds, res_min, res_max)
    try:
        aa_low, aa_high = lookup_gap_estimate(agg, gap_ang)
    except GapEstimateUnavailable as e:
        raise GapDataError(f"gap {label}: {e}") from e
    print(
        f"  Gap {label}: {gap_ang} \N{ANGSTROM SIGN}  ->  {aa_low}-{aa_high} residues",
        file=sys.stderr,
    )
    return aa_low, aa_high


def _chain_gap_estimates(chain, gaps, res_min, res_max, agg, thresholds):
    """
    Per-gap (aa_low, aa_high, prev_resnum, curr_resnum) for one chain's internal
    gaps, plus the raw Cα-Cα distances they were estimated from.
    """
    cid = chain.get_id()
    distances = ca_distances(chain)

    estimates, angstroms = [], []
    for (i_prev, i_curr, prev_resnum, curr_resnum) in gaps:
        gap_ang = math.floor(distances[i_prev][i_curr])
        label = f"{cid}{prev_resnum}-{curr_resnum}"
        aa_low, aa_high = _estimate(agg, gap_ang, label, res_min, res_max, thresholds)
        estimates.append((aa_low, aa_high, prev_resnum, curr_resnum))
        angstroms.append(gap_ang)
    return estimates, angstroms


def _gapped_chain_parts(chain, gap_estimates, term_lo, term_hi, sep):
    """One chain with internal gaps: [term/]fixed/gap/fixed[/.../term], joined by sep."""
    cid = chain.get_id()
    first_res, last_res = chain_span(chain)

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


def _chain_order_estimates(chains, chain_ids, res_min, res_max, agg, thresholds):
    """Per-gap (aa_low, aa_high) between consecutive chains, plus their Cα-Cα distances."""
    estimates, angstroms = [], []
    for i in range(len(chains) - 1):
        gap_ang = math.floor(interchain_gap_distance(chains[i], chains[i + 1]))
        label = f"{chain_ids[i]}(last)-{chain_ids[i + 1]}(first)"
        estimates.append(_estimate(agg, gap_ang, label, res_min, res_max, thresholds))
        angstroms.append(gap_ang)
    return estimates, angstroms


def _chain_order_parts(chains, gap_estimates, term_lo, term_hi, sep):
    """Chains merged into one designed chain: [term/]fixed/gap/fixed[/.../term], joined by sep."""
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


def n_motif_residues(chains):
    """Residues in the designed segments. Chains carried through unchanged don't count."""
    return sum(len(Selection.unfold_entities(c, "R")) for c in chains)


def build_contig(structure, res_min, res_max, gap_size_data, chain_order=None, sep=","):
    """
    Assemble the full contig body (without the surrounding "contigmap.contigs=[...]"
    or '"contig": "..."' wrapper): one designed/fixed segment per output
    chain, joined by a hard chain break ("/0 "). `sep` is the intra-chain
    separator: "," for the standard style, "/" for --rfd1.

    res_min/res_max are whole-design lengths, so the motif count and the terminal
    budget are pooled over EVERY designed segment and only then split back over
    the 2 * n_segments terminals. Budgeting per segment would credit each one
    with the residues of all the others -- with two gapped chains of 60 aa each,
    both would compute their terminals as if 60 aa of the requested length were
    free rather than 120. Chains carried through unchanged are not designed and
    count towards neither the motif total nor the budget.
    """
    all_chains = list(structure[0].get_list())
    print(f"  Chains found: {', '.join(c.get_id() for c in all_chains)}", file=sys.stderr)

    agg = aggregate_by_residue_range(gap_size_data, res_min, res_max)
    thresholds = gap_size_thresholds(gap_size_data)

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

        motif_residues = n_motif_residues(ordered_chains)
        check_motif_fits(motif_residues, res_min)

        gap_estimates, gap_angstroms = _chain_order_estimates(
            ordered_chains, chain_ids, res_min, res_max, agg, thresholds
        )
        print(file=sys.stderr)

        term_lo, term_hi = resolve_terminals(
            motif_residues,
            sum(lo for lo, _ in gap_estimates),
            sum(hi for _, hi in gap_estimates),
            res_min, res_max, n_segments=1,
            gap_size_data=gap_size_data, gap_angstroms=gap_angstroms,
        )
        motif_segment = _chain_order_parts(
            ordered_chains, gap_estimates, term_lo, term_hi, sep
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

    designed = [chain for chain in all_chains if chain.get_id() in gaps_by_id]
    motif_residues = n_motif_residues(designed)
    check_motif_fits(motif_residues, res_min)

    estimates_by_id, gap_angstroms = {}, []
    for chain in designed:
        cid = chain.get_id()
        estimates, angstroms = _chain_gap_estimates(
            chain, gaps_by_id[cid], res_min, res_max, agg, thresholds
        )
        estimates_by_id[cid] = estimates
        gap_angstroms.extend(angstroms)
    print(file=sys.stderr)

    all_estimates = [e for estimates in estimates_by_id.values() for e in estimates]
    term_lo, term_hi = resolve_terminals(
        motif_residues,
        sum(lo for lo, hi, *_ in all_estimates),
        sum(hi for lo, hi, *_ in all_estimates),
        res_min, res_max, n_segments=len(designed),
        gap_size_data=gap_size_data, gap_angstroms=gap_angstroms,
    )

    groups = []
    for chain in all_chains:
        cid = chain.get_id()
        if cid in estimates_by_id:
            groups.append(
                _gapped_chain_parts(chain, estimates_by_id[cid], term_lo, term_hi, sep)
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
    parser.add_argument(
        "res_max", type=int, nargs="?", default=None,
        help=(
            "Maximum total protein length (10-499). If omitted, res_min is treated "
            "as a single fixed target length and \"length\" is printed as that single "
            "value instead of a min-max range."
        ),
    )
    parser.add_argument(
        "--pickle-file",
        default=DEFAULT_CHECKPOINT,
        help=(
            "Checkpoint dataset to use: the bundled 'gyr' variant (default) "
            "or a path to an external .pkl file."
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

    fixed_length = args.res_max is None
    res_max = args.res_min if fixed_length else args.res_max

    if not (10 <= args.res_min <= 499 and 10 <= res_max <= 499):
        print("Error: res_min and res_max must both be in [10, 499].", file=sys.stderr)
        sys.exit(1)

    if args.res_min > res_max:
        print("Error: res_min must be <= res_max.", file=sys.stderr)
        sys.exit(1)

    try:
        gap_size_data = load_pickle(args.pickle_file)
    except (OSError, ValueError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    structure = load_pdb(args.pdb_file)

    try:
        body = build_contig(
            structure, args.res_min, res_max, gap_size_data,
            chain_order=args.chain_order, sep="/" if args.rfd1 else ",",
        )
    except AutoContigmapError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

    length_range = str(args.res_min) if fixed_length else f"{args.res_min}-{res_max}"
    if args.rfd1:
        print(f"contigmap.contigs=[{body}]")
        print(f'contigmap.length="{length_range}"')
    else:
        print(f'"contig": "{body}"')
        print(f'"length": "{length_range}"')


if __name__ == "__main__":
    main()
