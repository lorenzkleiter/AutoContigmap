"""
core.py — motif gap-size estimation from Cα-Cα distances.

Gap sizes are estimated by looking the target Cα-Cα distance up in one of the
bundled residue-count-vs-distance checkpoints (data/ — built from real PDB
statistics), aggregated over a caller-supplied total-protein-length range.

Library entry points: estimate_gap_fill() for a single distance,
gap_size_percentile_threshold() for the extreme-distance warning threshold.
Contig assembly lives in cli.py.
"""

import math
import pickle
from importlib import resources

import numpy as np

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection
from Bio.PDB.internal_coords import AtomKey
from Bio.PDB.ic_rebuild import structure_rebuild_test

DEFAULT_CHECKPOINT = "results_checkpoint_gyr"
MIN_CHECKPOINT_RESIDUES = 10
MAX_CHECKPOINT_RESIDUES = 499

# A bonded Cα-Cα pair sits at ~3.8 Å; beyond this the backbone is broken.
CA_GAP_DISTANCE_THRESHOLD = 4.0

GAP_SIZE_WARN_PERCENTILE = 0.95

# Past this percentile of the pair-distance distribution an estimate rests on
# almost no data, so it is a hard error rather than a warning.
GAP_SIZE_ERROR_PERCENTILE = 0.99

# Per-length thresholds are estimated from one chain length's own pairs, so they
# are noisy: on a real checkpoint 215 of 489 consecutive steps run downhill, with
# single-length dips of up to 37 Å. Longer chains can only reach further, so the
# true threshold is non-decreasing in length; a rolling median followed by a
# running maximum removes the dips without inventing structure. Without it, one
# noisy length in the middle of a window rejects a gap that every other length
# accepts.
THRESHOLD_SMOOTHING_WINDOW = 11

QUANTILE_KEYS = ("Q0.25", "Q0.4", "Q0.5", "Q0.6", "Q0.75")


class GapEstimateUnavailable(Exception):
    """Raised when a gap distance has no usable estimate in the checkpoint."""


def load_pdb(pdb):
    parser = PDBParser(PERMISSIVE=0, QUIET=True)
    return parser.get_structure("motif", pdb)


def load_pickle(name: str = DEFAULT_CHECKPOINT):
    """
    Load a residue-count-vs-Cα-distance checkpoint dataset.

    `name` is either one of the bundled variants (results_checkpoint_gyr,
    results_checkpoint_gyr_ss_2, results_checkpoint_standard,
    results_checkpoint_surface — with or without the "results_checkpoint_"
    prefix) or a filesystem path to an external .pkl file.
    """
    bundled = resources.files("autocontigmap.data") / f"{name}.pkl"
    if bundled.is_file():
        with bundled.open("rb") as f:
            return pickle.load(f)

    short = resources.files("autocontigmap.data") / f"results_checkpoint_{name}.pkl"
    if short.is_file():
        with short.open("rb") as f:
            return pickle.load(f)

    with open(name, "rb") as f:
        return pickle.load(f)


# ---------------------------------------------------------------------------
# Aggregation over the requested design-length window
# ---------------------------------------------------------------------------

def aggregate_by_residue_range(Results, res_min, res_max):
    """
    Aggregate the checkpoint matrices over a range of chain lengths.

    Quantiles are averaged WITHOUT weights: every admissible design length in
    [res_min, res_max] counts equally. Weighting by chain counts would import
    the length distribution of the PDB, which has nothing to do with what the
    user asked for; weighting by pair counts would upweight long chains by a
    factor growing with n(n-1)/2.

    'counts' keeps chain-count weighting, because it feeds the gap-size
    percentile threshold, where pooled pairs per chain over the actual
    population is the quantity of interest.

    Cells with no data are NaN in the checkpoint and are simply skipped: the
    mean is over however many chain lengths carry data for that distance bin,
    and the bin is NaN only when none do. No coverage threshold is applied here
    — distances where coverage is poor are the distances the per-length
    percentile check (check_gap_against_lengths) has already rejected, so those
    rows are never read. '<key>_coverage' is still returned, purely as a
    diagnostic.

    Parameters
    ----------
    Results  : dict — the full checkpoint dictionary
    res_min  : int  — minimum chain length, inclusive (10-499)
    res_max  : int  — maximum chain length, inclusive (10-499)

    Returns
    -------
    agg : dict with matrix-valued keys reduced to 1-D arrays over the distance
          axis, quantile keys gaining a '<key>_coverage' companion, and
          'no_pairs_aggregated' reduced to the total chain count in the window.
          Non-matrix keys pass through unchanged.
    """
    if res_min > res_max:
        raise ValueError(f"res_min ({res_min}) must not exceed res_max ({res_max}).")
    if not (MIN_CHECKPOINT_RESIDUES <= res_min <= MAX_CHECKPOINT_RESIDUES) or not (
        MIN_CHECKPOINT_RESIDUES <= res_max <= MAX_CHECKPOINT_RESIDUES
    ):
        raise ValueError(
            f"Residue range [{res_min}, {res_max}] is outside "
            f"[{MIN_CHECKPOINT_RESIDUES}, {MAX_CHECKPOINT_RESIDUES}]."
        )
    shape = Results["Q0.5"].shape
    n_dist = shape[0]
    col_min = res_min - MIN_CHECKPOINT_RESIDUES
    col_max = res_max - MIN_CHECKPOINT_RESIDUES + 1
    n_window = col_max - col_min

    weights_full = Results["no_pairs_aggregated"]  # chains per length, broadcast down rows

    agg = {}
    for key, val in Results.items():
        if not (isinstance(val, np.ndarray) and val.shape == shape):
            agg[key] = val
            continue

        block = val[:, col_min:col_max]

        if key in QUANTILE_KEYS:
            valid = ~np.isnan(block)
            n_valid = valid.sum(axis=1)
            coverage = n_valid / n_window
            agg[key] = np.divide(
                np.where(valid, block, 0.0).sum(axis=1),
                n_valid,
                out=np.full(n_dist, np.nan),
                where=n_valid > 0,
            )
            agg[key + "_coverage"] = coverage

        elif key == "counts":
            w = np.where(np.isnan(block), 0.0, weights_full[:, col_min:col_max])
            wsum = w.sum(axis=1)
            agg[key] = np.divide(
                (w * np.nan_to_num(block)).sum(axis=1),
                wsum,
                out=np.full(n_dist, np.nan),
                where=wsum > 0,
            )

        elif key == "no_pairs_aggregated":
            # total chains in the window; the stored matrix is constant down rows
            agg[key] = float(block[0, :].sum())

        elif key == "cumsum":
            continue  # recomputed on demand from the aggregated counts

        else:
            agg[key] = np.nanmean(block, axis=1)

    return agg


def _clamp_residue_range(res_min, res_max):
    """Clamp to the checkpoint's covered [10, 499] range; raise if there's no overlap at all."""
    if res_max < MIN_CHECKPOINT_RESIDUES or res_min > MAX_CHECKPOINT_RESIDUES:
        raise ValueError(
            f"Residue range [{res_min}, {res_max}] doesn't overlap the checkpoint's "
            f"covered range [{MIN_CHECKPOINT_RESIDUES}, {MAX_CHECKPOINT_RESIDUES}]."
        )
    return max(MIN_CHECKPOINT_RESIDUES, int(res_min)), min(MAX_CHECKPOINT_RESIDUES, int(res_max))


# ---------------------------------------------------------------------------
# Distance-bin lookup
# ---------------------------------------------------------------------------

def distance_bin(gap_ang, n_dist):
    """
    Row index for a Cα-Cα distance. Distance bins start at 1-1.999 Å, so row r
    covers [r + 1, r + 2) Å and a distance d maps to floor(d) - 1.
    """
    row = int(math.floor(gap_ang)) - 1
    if row < 0 or row >= n_dist:
        raise GapEstimateUnavailable(
            f"Cα-Cα distance {gap_ang:.1f} Å is outside the checkpoint's binned "
            f"range [1, {n_dist + 1}) Å."
        )
    return row


def lookup_gap_estimate(agg, gap_ang):
    """
    (aa_low, aa_high) residue-count estimate for a gap of `gap_ang` Å, read off
    the aggregated Q0.4 / Q0.6 curves.

    Raises GapEstimateUnavailable only when no estimate exists at all: the
    distance is outside the binned range, or no chain length in the window has
    residue pairs that far apart. In practice the per-length percentile check
    (check_gap_against_lengths) rejects such distances first, so this is a
    backstop rather than the usual path.
    """
    row = distance_bin(gap_ang, len(agg["Q0.4"]))
    q_low, q_high = agg["Q0.4"][row], agg["Q0.6"][row]

    if np.isnan(q_low) or np.isnan(q_high):
        raise GapEstimateUnavailable(
            f"No gap estimate for a {math.floor(gap_ang):.0f} Å Cα-Cα distance in this "
            f"design-length window: no chain length in {len(agg['Q0.4'])}-bin data has "
            f"residue pairs that far apart. Raise or widen the length range."
        )

    # The checkpoint quantiles are over sequence SEPARATIONS (r_j - r_i), while a
    # contig declares the number of residues INSERTED between the flanking
    # residues, hence the -1 on both bounds. Rounding is outward so the interval
    # widens rather than narrows: ceil on the lower bound, floor on the upper.
    aa_low = max(0, int(math.ceil(q_low)) - 1)
    aa_high = max(0, int(math.floor(q_high)) - 1)
    # Both quantiles can land inside the same integer, which would inverse the
    # interval; a gap can never require fewer residues than its own lower bound.
    return aa_low, max(aa_low, aa_high)


def length_distance_thresholds(Results, percentile):
    """
    Per chain length, the Cα-Cα distance (Å) at `percentile` of THAT length's own
    pair-distance distribution. Returns a 1-D array over chain lengths
    10..MAX_CHECKPOINT_RESIDUES, NaN where a length has no data.

    This is deliberately computed before any aggregation over the requested
    design-length window. The pooled distribution is dominated by the longest
    chains in the window, because pairs per chain grow as n(n-1)/2, so a single
    aggregated threshold silently under-warns for the short end of the range.
    Per length, a gap can be flagged as implausible for 90-residue designs and
    perfectly ordinary for 110-residue ones, which is both the truth and the
    actionable form: it says what res_min would have to be.

    The threshold is the upper edge of the first distance bin whose cumulative
    share exceeds `percentile`, matching the [r + 1, r + 2) Å binning, then
    smoothed and made non-decreasing (see THRESHOLD_SMOOTHING_WINDOW).
    """
    counts = np.nan_to_num(np.asarray(Results["counts"]))
    total = counts.sum(axis=0)
    with np.errstate(invalid="ignore", divide="ignore"):
        cdf = np.cumsum(counts, axis=0) / np.where(total > 0, total, np.nan)
    over = cdf > percentile
    any_over = over.any(axis=0)
    thresholds = np.where(any_over, np.argmax(over, axis=0) + 2.0, np.nan)
    thresholds = np.where(total > 0, thresholds, np.nan)
    return _smooth_monotone(thresholds, THRESHOLD_SMOOTHING_WINDOW)


def _smooth_monotone(values, window):
    """Rolling median, then enforce non-decreasing. NaN-tolerant; all-NaN stays NaN."""
    if np.isnan(values).all():
        return values
    pad = window // 2
    padded = np.pad(values, (pad, pad), mode="edge")
    with np.errstate(invalid="ignore"):
        med = np.array([np.nanmedian(padded[i:i + window]) for i in range(len(values))])
    filled = np.where(np.isnan(med), -np.inf, med)
    out = np.maximum.accumulate(filled)
    return np.where(np.isinf(out), np.nan, out)


class GapSizeVerdict:
    """
    Outcome of checking one gap distance against every length in the window.

    level        : "ok" | "warn" | "error"
    bad_lengths  : the chain lengths for which the distance is beyond the
                   relevant percentile, as a (first, last) pair, or None
    n_bad        : how many lengths those are
    n_total      : lengths in the requested window
    min_viable   : smallest length from which onwards EVERY length in the
                   window clears the threshold, or None if no length does
    threshold    : the percentile threshold at min_viable, for the message
    """

    def __init__(self, level, bad_lengths, n_bad, n_total, min_viable, percentile):
        self.level = level
        self.bad_lengths = bad_lengths
        self.n_bad = n_bad
        self.n_total = n_total
        self.min_viable = min_viable
        self.percentile = percentile


def _verdict(gap_ang, lengths, thresholds, percentile, level):
    # NaN threshold means the length has no data at all, which counts as "does
    # not support this distance" rather than silently passing.
    bad = ~(gap_ang <= thresholds)
    if not bad.any():
        return None
    good = ~bad
    # smallest length from which onwards everything clears, so a noisy
    # threshold in the middle of the range cannot suggest an unusable res_min
    ok_from = None
    for i in range(len(lengths)):
        if good[i:].all():
            ok_from = int(lengths[i])
            break
    bad_ls = lengths[bad]
    return GapSizeVerdict(
        level, (int(bad_ls[0]), int(bad_ls[-1])), int(bad.sum()),
        len(lengths), ok_from, percentile,
    )


def check_gap_against_lengths(gap_ang, res_min, res_max, thr_warn, thr_error):
    """
    Check one gap distance against every chain length in [res_min, res_max].

    Returns a GapSizeVerdict with level "error" if the distance is beyond
    GAP_SIZE_ERROR_PERCENTILE for any length in the window, "warn" if beyond
    GAP_SIZE_WARN_PERCENTILE for any length, otherwise None.
    """
    lengths = np.arange(res_min, res_max + 1)
    cols = lengths - MIN_CHECKPOINT_RESIDUES
    err = _verdict(gap_ang, lengths, thr_error[cols], GAP_SIZE_ERROR_PERCENTILE, "error")
    if err is not None:
        return err
    return _verdict(gap_ang, lengths, thr_warn[cols], GAP_SIZE_WARN_PERCENTILE, "warn")


def estimate_gap_fill(gap_ang, res_min, res_max, checkpoint=DEFAULT_CHECKPOINT, gap_size_data=None):
    """
    (aa_low, aa_high) residue-count estimate for a gap of gap_ang Angstrom,
    aggregated over the res_min-res_max total-protein-length range (clamped
    to the checkpoint's covered [10, 499] range). Pass a pre-loaded
    gap_size_data (from load_pickle()) to avoid re-reading the checkpoint
    file for every call.
    """
    if gap_size_data is None:
        gap_size_data = load_pickle(checkpoint)
    res_min, res_max = _clamp_residue_range(res_min, res_max)
    agg = aggregate_by_residue_range(gap_size_data, res_min, res_max)
    return lookup_gap_estimate(agg, gap_ang)


# ---------------------------------------------------------------------------
# Terminal budget
# ---------------------------------------------------------------------------

def pooled_terminal_budget(motif_residues, gap_lo_sum, gap_hi_sum, res_min, res_max):
    """
    Permissive pooled terminal budget and its feasibility status.

        aa_term_low  = res_min - motif - gap_hi_sum
        aa_term_high = res_max - motif - gap_lo_sum

    These are the projection of the feasible set onto the terminal axis: every
    value in between is realisable by some admissible gap draw, and nothing
    outside is. They do NOT by themselves pin [res_min, res_max] — the emitted
    contig also needs the length flag.

    Returns (aa_term_low, aa_term_high, status), status one of:

      "ok"              both bounds non-negative.
      "clamped"         aa_term_low < 0 but res_min is still reachable with the
                        gaps near their lower estimates; aa_term_low returns 0.
      "unreachable_min" res_min < motif + gap_lo_sum, so even the shortest
                        buildable object overshoots res_min (equivalently
                        aa_term_low + sigma_gap < 0). The realised minimum is
                        motif + gap_lo_sum.
      "infeasible"      aa_term_high < 0, i.e. res_max < motif + gap_lo_sum.
                        No admissible configuration exists at all.
    """
    aa_term_low = res_min - motif_residues - gap_hi_sum
    aa_term_high = res_max - motif_residues - gap_lo_sum
    floor_length = motif_residues + gap_lo_sum  # shortest buildable object

    if aa_term_high < 0:
        return max(0, aa_term_low), aa_term_high, "infeasible"
    if aa_term_low < 0:
        # aa_term_low + sigma_gap == res_min - motif - gap_lo_sum
        status = "clamped" if res_min >= floor_length else "unreachable_min"
        return 0, aa_term_high, status
    return aa_term_low, aa_term_high, "ok"


def split_terminal_budget(aa_term_low, aa_term_high, n_segments=1):
    """
    Split a pooled budget over the terminals actually emitted.

    `n_segments` is the number of designed segments that each receive an N- and
    a C-terminal range, so the budget divides over 2 * n_segments terminals.
    Both bounds round UP: with an even budget the ceiling is inert, and with an
    odd one it inflates the declared pooled range by one residue rather than
    deciding a priori which terminal carries the extra residue. That surplus is
    not realised in practice, because the length flag pins [res_min, res_max]
    regardless of what the declared ranges permit.
    """
    n_term = 2 * max(1, n_segments)
    lo = max(0, math.ceil(aa_term_low / n_term))
    hi = max(0, math.ceil(aa_term_high / n_term))
    return lo, max(lo, hi)


# ---------------------------------------------------------------------------
# Cα distances and gap detection
# ---------------------------------------------------------------------------

def ca_distances(chain):
    structure_rebuild_test(chain)
    atm_name_ndx = AtomKey.fields.atm
    aa_i = chain.internal_coord.atomArrayIndex
    ca_select = [aa_i.get(k) for k in aa_i.keys() if k.akl[atm_name_ndx] == "CA"]
    return chain.internal_coord.distance_plot(ca_select)


def find_gaps_in_chain(chain, distance_threshold=CA_GAP_DISTANCE_THRESHOLD):
    """
    [(i_prev, i_curr, prev_resnum, curr_resnum), ...] for every internal gap.

    A gap is flagged on a residue-numbering jump OR on a Cα-Cα distance beyond
    `distance_threshold`, since some structures break the backbone while
    keeping sequential numbering (renumbered chains) and others jump in
    numbering without an actual break (insertion codes, renumbering schemes).
    """
    res_list = Selection.unfold_entities(chain, "R")
    gaps = []
    for i in range(1, len(res_list)):
        prev_res, curr_res = res_list[i - 1], res_list[i]
        prev_num = prev_res.get_id()[1]
        curr_num = curr_res.get_id()[1]
        numbering_gap = curr_num - prev_num > 1
        ca_gap = (
            "CA" not in prev_res
            or "CA" not in curr_res
            or np.linalg.norm(prev_res["CA"].coord - curr_res["CA"].coord) > distance_threshold
        )
        if numbering_gap or ca_gap:
            gaps.append((i - 1, i, prev_num, curr_num))
    return gaps
