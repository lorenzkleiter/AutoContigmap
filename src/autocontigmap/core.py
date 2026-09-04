"""
core.py — motif gap-size estimation from Cα-Cα distances, and the
higher-level chain_analyser()/contigmap() convenience API for building an
RFdiffusion contig string directly from a motif PDB.

Gap sizes are estimated by looking the target Cα-Cα distance up in one of
the bundled residue-count-vs-distance checkpoints (data/ — built from real
PDB statistics), aggregated over a caller-supplied total-protein-length
range.
"""

import math
import pickle
from collections import defaultdict
from importlib import resources

import numpy as np

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection
from Bio.PDB.internal_coords import AtomKey
from Bio.PDB.ic_rebuild import structure_rebuild_test

DEFAULT_CHECKPOINT = "results_checkpoint_gyr"
MIN_CHECKPOINT_RESIDUES = 10
MAX_CHECKPOINT_RESIDUES = 499


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


def aggregate_by_residue_range(Results, res_min, res_max):
    """
    Aggregate Results matrices over a range of residue lengths using
    weighted averaging, where weights are Results['no_pairs_aggregated'].

    Parameters
    ----------
    Results  : dict  — the full results dictionary (248 x 490 matrices)
    res_min  : int   — minimum residue length, inclusive (10-499)
    res_max  : int   — maximum residue length, inclusive (10-499)

    Returns
    -------
    agg : dict with the same keys, aggregated over the residue axis.
          Matrix-valued keys become 1-D arrays of length 248.
          Scalar/non-matrix keys are passed through unchanged.
    """
    col_min = res_min - 10
    col_max = res_max - 10 + 1

    if col_min < 0 or col_max > 490:
        raise ValueError(f"Residue range [{res_min}, {res_max}] is outside [10, 499].")

    quantile_keys = {"Q0.25", "Q0.4", "Q0.5", "Q0.6", "Q0.75"}
    weighted_keys = {"counts", "cumsum", "no_pairs_aggregated"}

    agg = {}
    for key, val in Results.items():
        if isinstance(val, np.ndarray) and val.shape == (248, 490):
            slice_ = val[:, col_min:col_max].copy()

            if key in quantile_keys:
                w = Results["no_pairs_aggregated"][:, col_min:col_max].copy()
                nan_mask = np.isnan(slice_)
                slice_[nan_mask] = 0.0
                w[nan_mask] = 0.0

                weight_sum = w.sum(axis=1)
                safe_sum = np.where(weight_sum > 0, weight_sum, 1)
                agg[key] = (w * slice_).sum(axis=1) / safe_sum

            elif key in weighted_keys:
                w = Results["no_pairs_aggregated"][:, col_min:col_max]
                weight_sum = w.sum(axis=1)
                safe_sum = np.where(weight_sum > 0, weight_sum, 1)
                agg[key] = (w * slice_).sum(axis=1) / safe_sum

            else:
                agg[key] = slice_.mean(axis=1)
        else:
            agg[key] = val

    return agg


def _clamp_residue_range(res_min, res_max):
    """Clamp to the checkpoint's covered [10, 499] range; raise if there's no overlap at all."""
    if res_max < MIN_CHECKPOINT_RESIDUES or res_min > MAX_CHECKPOINT_RESIDUES:
        raise ValueError(
            f"Residue range [{res_min}, {res_max}] doesn't overlap the checkpoint's "
            f"covered range [{MIN_CHECKPOINT_RESIDUES}, {MAX_CHECKPOINT_RESIDUES}]."
        )
    return max(MIN_CHECKPOINT_RESIDUES, int(res_min)), min(MAX_CHECKPOINT_RESIDUES, int(res_max))


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
    gap_ang = int(math.floor(gap_ang))
    aa_low = int(np.floor(agg["Q0.4"][gap_ang - 1]))
    aa_high = int(np.floor(agg["Q0.6"][gap_ang - 1]))
    return aa_low, aa_high


def gap_size_percentile_threshold(res_min, res_max, percentile=0.95, checkpoint=DEFAULT_CHECKPOINT, gap_size_data=None):
    """
    Gap size (Å) at `percentile` of the gap-size distribution seen in the
    checkpoint data for the res_min-res_max range (clamped to [10, 499]),
    or None if it can't be determined (degenerate distribution).
    """
    if gap_size_data is None:
        gap_size_data = load_pickle(checkpoint)
    res_min, res_max = _clamp_residue_range(res_min, res_max)
    agg = aggregate_by_residue_range(gap_size_data, res_min, res_max)
    a = np.cumsum(agg["counts"])
    norm = (a - a.min()) / (a.max() - a.min())
    idx = np.where(norm > percentile)[0]
    return int(idx[0]) + 1 if len(idx) else None


def ca_distances(chain):
    structure_rebuild_test(chain)
    atm_name_ndx = AtomKey.fields.atm
    aa_i = chain.internal_coord.atomArrayIndex
    ca_select = [aa_i.get(k) for k in aa_i.keys() if k.akl[atm_name_ndx] == "CA"]
    return chain.internal_coord.distance_plot(ca_select)


def sizetest(agg, ca_dist):
    a = np.cumsum(agg["counts"])
    norm = (a - a.min()) / (a.max() - a.min())
    itemindex = np.where(norm > 0.90)
    return np.max(ca_dist) > itemindex[0][0]


def gap_estimation(chain, gap, res, gap_size_data=None):
    """
    chain : Bio.PDB chain
    gap   : [i_prev, i_curr] — residue-list indices (not residue numbers)
    res   : [res_min, res_max]
    gap_size_data : pre-loaded checkpoint dict, or None to load the default

    Returns (aa_low, aa_high, gap_angstrom).
    """
    if gap_size_data is None:
        gap_size_data = load_pickle()

    ca_dist = ca_distances(chain)
    # gap size is the floor of ca_dist at the gap (distance bins start at 1-1.999 Angstrom)
    gap_size = math.floor(ca_dist[gap[0]][gap[1]])
    agg = aggregate_by_residue_range(gap_size_data, res[0], res[1])

    if sizetest(agg, ca_dist):
        print(f"Warning your maximum pair distance in your motif is {np.max(ca_dist)}")
        print(f"Such a distance is very unlikely in a protein of size {res[0]} - {res[1]}")
        print("Residue pair distance data might be unreliable")

    # gap sizes are stored one lower and residues 10 higher in the data
    aa_low = int(np.floor(agg["Q0.4"][gap_size - 1]))
    aa_high = int(np.floor(agg["Q0.6"][gap_size - 1]))
    return aa_low, aa_high, gap_size


def chain_analyser(pdb, res, gap_size_data=None):
    """
    pdb : path to pdb
    res : [res_min, res_max]
    """
    if gap_size_data is None:
        gap_size_data = load_pickle()

    structure = load_pdb(pdb)
    chain_information = defaultdict()

    for j, chain in enumerate(structure[0].get_list()):
        res_list = Selection.unfold_entities(chain, "R")
        temp_dict = {
            "id": chain.get_id(),
            "first_res": res_list[0].get_id()[1],
            "last_res": res_list[len(res_list) - 1].get_id()[1],
            "gaps": [],
            "gap_sizes": [],
            "max_size": len(res_list),
            "motif_size": len(res_list),
            "min_size": len(res_list),
            "res": res,
        }
        for i in range(1, len(res_list)):
            prev_res = res_list[i - 1]
            curr_res = res_list[i]
            prev_num = prev_res.get_id()[1]
            curr_num = curr_res.get_id()[1]

            if curr_num - prev_num > 1:
                temp_dict["gaps"].append((prev_num, curr_num))
                gap_stats = gap_estimation(chain, [i - 1, i], res, gap_size_data)
                temp_dict["gap_sizes"].append(gap_stats)
                temp_dict["max_size"] += gap_stats[1]
                temp_dict["min_size"] += gap_stats[0]
        chain_information[j] = temp_dict

    return chain_information


def contigmap(pdb, res, gap_size_data=None):
    """
    pdb : file path to pdb
    res : [res_min, res_max]

    Returns an RFdiffusion contig string, or "" if the requested length
    range can't accommodate the estimated gaps.
    """
    too_small = False
    chains_info = chain_analyser(pdb, res, gap_size_data)
    string = "["
    for i, chain in enumerate(chains_info):
        if chains_info[i]["max_size"] > chains_info[i]["res"][1] and chains_info[i]["gaps"] != []:
            print("The upper limit of the amount of residues necessary for the gaps + the amount in the motif is bigger than the upper range given")
            print("You need a higher upper limit of residues for", chains_info[i]["id"])
            print("The gaps in your motif are of size:", chains_info[i]["gap_sizes"])
            too_small = True
            break

        elif chains_info[i]["min_size"] > chains_info[i]["res"][0] and chains_info[i]["gaps"] != []:
            print("The lower limit of the amount of residues necessary for the gaps + the amount in the motif is bigger than the lower range given")
            print("You need a smaller lower limit of residues for", chains_info[i]["id"])
            print("The gaps in your motif are of size:", chains_info[i]["gap_sizes"])
            too_small = True
            break

        if chains_info[i]["gaps"] != [] or i == 0:
            left_over_high = math.floor((chains_info[i]["res"][1] - chains_info[i]["max_size"]) / 2)
            string += f"0-{left_over_high}/"
        string += chains_info[i]["id"]
        string += str(chains_info[i]["first_res"])
        string += "-"

        if chains_info[i]["gaps"] != []:
            for j, gap in enumerate(chains_info[i]["gaps"]):
                string += str(gap[0])
                string += f"/{chains_info[i]['gap_sizes'][j][0]}-{chains_info[i]['gap_sizes'][j][1]}/"
                string += chains_info[i]["id"]
                string += str(gap[1])
                string += "-"

        string += str(chains_info[i]["last_res"])

        if chains_info[i]["gaps"] != [] or i == 0:
            string += f"/0-{left_over_high}"

        if i < len(chains_info) - 1:
            string += "/0 "
        print(f"Its lower size range is {chains_info[i]['min_size']} and its maximum {chains_info[i]['max_size']} before filling")
    string += "]"

    if too_small:
        return ""
    return string
