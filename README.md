# AutoContigmap

Estimate RFdiffusion contig gap sizes for a motif PDB, from PDB
residue-count-vs-Cα-distance statistics, instead of guessing a gap/length
range by hand.

Given a motif PDB with one or more chain breaks (or several single-segment
chains meant to be scaffolded together) and a target total protein-length
range, it looks up each gap's Cα-Cα distance against a bundled statistics
checkpoint and reports a plausible residue-count range per gap, plus N-/C-
terminal tail budgets as an ready RFD3 or RFD contigmap string.

## Install

```bash
pip install "git+https://github.com/lorenzkleiter/AutoContigmap.git"
```

## CLI usage

```bash
autocontigmap motif.pdb <res_min> [<res_max>]
```

If `res_max` is omitted, `res_min` is treated as a single fixed target
length (internally `res_max = res_min`) and `length` is printed as that
single value (`"408"`) rather than a `min-max` range.

Note: motif.pdb should only include the motif residues. Gaps are detected from that automatically. There is no way to input the whole protein and define the motif afterwards.

Output on stdout, with fixed motif/target chain spans included alongside
the estimated gap-fill ranges:

```
"contig": "8-15,A18-25,16-30,A47-54,16-30,A92-99"
"length": "150-200"
```

`length` is just `res_min-res_max` (or the single `res_min`) echoed back,
for pinning the overall
design length directly (the default terminal budget below otherwise only
pins a `0-upper_bound` range, not an exact total).

Pass `--rfd1` to instead print the equivalent old-style RFdiffusion(1)
contig, slash-separated:

```
contigmap.contigs=[8-15/A18-25/16-30/A47-54/16-30/A92-99]
contigmap.length="150-200"
```

A terminal segment that computes to exactly `0-0` is omitted from the
contig entirely rather than printed as a literal `0`.

If a gap's estimated Cα-Cα distance is beyond the 95th percentile of gap
sizes seen in the checkpoint data for the requested `res_min`-`res_max`
range, a warning is printed to stderr — the estimate is based on thin data
out there, and a larger design is probably necessary.

### Diagnostics (stderr)

Every run starts by reporting what it found, before anything else:

```
  Chains found: A, B
  Motif identified automatically (has internal gap(s)): A
```

(or, with `--chain-order`, `Motif identified via --chain-order: B, A (merged into one chain, gaps measured between consecutive chains)`.)

Then one `Gap ...: X Å -> lo-hi residues` line per gap, each immediately
followed by its own 95th-percentile warning if it applies; a blank line;
then the length-range warning below if it applies. All diagnostic lines use
the same two-space indent and `Warning: ...` phrasing, and both warnings
print if both apply. Hard errors are checked before any of this is printed,
so an errored run never prints a warning it's about to make moot.

## Length-range sanity checks

- **Hard error** if `res_min` is smaller than the motif's own residue count
  (excluding gaps) — the requested range can't even fit the fixed motif
  residues, let alone the gaps. Checked before any gap processing, so
  nothing else prints beforehand:

  ```
  Error: res_min (20) is smaller than the motif's own residue count (32) -- the requested length range can't even fit the fixed motif residues, let alone the gaps.
  ```

- **Hard error** if `res_min` is below `motif + gap_lo_sum` — the lower end
  of the requested range can't be reached at any terminal length:

  ```
  Error: res_min (181) is below the shortest buildable design (120 aa motif + 100-133 aa of gaps) -- try res_min >= 225, or >= 268 (recommended).
  ```

  A `res_max` below the same floor reports the same error. It is the same
  problem: `res_min <= res_max`, so a `res_max` under the floor puts
  `res_min` under it too, and raising `res_min` above it makes `res_max`
  fine automatically — so only `res_min` is ever suggested.

  The two values come from a scan that **re-estimates the gaps at every
  candidate**, against the `res_max` you gave. Sequence separation grows
  with chain length, so the `motif + gap_lo_sum` floor computed for the
  window that just failed is not itself a usable `res_min` — raising
  `res_min` to it re-estimates the gaps upward and moves the floor again.
  The first value is the smallest `res_min` that builds at all; the second,
  the recommended one, is the smallest at which the *high* end of every gap
  estimate still fits, so it also avoids the pinned-terminal warning below.
  If no `res_min` works at that `res_max`, the message says so instead:

  ```
  Error: res_min (130) is below the shortest buildable design (120 aa motif + 55-75 aa of gaps) -- res_min and res_max both have to be increased.
  ```

- **Warning** (not fatal) if the lower terminal budget is negative but
  every requested length is still reachable with the gaps near their lower
  estimates — the budget is pinned to 0 and the contig is still printed:

  ```
  Warning: the lower terminal budget is negative and was pinned to 0; res_min (230) is still reachable, but only when the gaps sample near their lower estimates (102 aa total).
  ```

Options:

- `--pickle-file NAME_OR_PATH` — which statistics checkpoint to use: the
  bundled `gyr` variant (default), or a path to an external `.pkl` file.
- `--rfd1` — output the old-style `contigmap.contigs=[...]`/`contigmap.length=...`
  form instead of the standard `"contig": ...`/`"length": ...` one.
- `--chain-order B,A` — for motifs whose segments are split across separate
  PDB chains (one chain per segment) instead of one chain with internal
  chain breaks. Gaps are then measured between the last residue of each
  chain and the first residue of the next.

## Chain selection

Any chain that contains at least one internal gap (chain break) is treated
as a motif chain and gap-filled. If no chain has a gap and `--chain-order`
isn't given, this is an error — there's nothing to scaffold:

```
Error: Not motif chain identified, motif chains have to have gaps. If the
motif is defined over multiple chains use --chain-order
```

**Default mode (no `--chain-order`)**: every chain with an internal gap is
gap-filled and becomes its own designed segment. Chains without a gap are
carried through unchanged as fixed spans. All of these segments are joined
with a hard chain break (`/0 `) in the PDB's chain order — each ends up as
its own separate output chain.

The terminal budget is **shared, not per segment**. `res_min`/`res_max` are
whole-design lengths, so the motif residues and gap estimates are pooled
over every designed segment and that single budget is split across all
`2 * n_segments` terminals. Chains carried through unchanged are not
designed and count towards neither.

**`--chain-order` mode**: the listed chains are instead scaffolded into a
single continuous designed chain (no break between them, one shared
terminal budget), with gaps measured between consecutive chains. Any other
chains present in the PDB are still carried through unchanged as fixed
spans, each behind its own chain break.

## Python API

The contig assembly itself lives in the CLI (`autocontigmap.cli`); the
package exports the pieces it is built from. For callers that already have
their own Cα-Cα distance measurement and just want the lookup:

```python
from autocontigmap import estimate_gap_fill, load_pickle

data = load_pickle()  # load once, reuse across calls
aa_low, aa_high = estimate_gap_fill(gap_ang=41, res_min=154, res_max=174, gap_size_data=data)
```

`estimate_gap_fill()` clamps `res_min`/`res_max` to the checkpoint's
covered `[10, 499]` range and raises `ValueError` if the requested range
doesn't overlap it at all.

To reproduce the CLI's plausibility check, `length_distance_thresholds()`
gives the per-chain-length Cα-Cα distance threshold at a percentile, and
`check_gap_against_lengths()` turns a distance into a `GapSizeVerdict`
(`level`, the lengths responsible, and the shortest design length that
would support the gap):

```python
from autocontigmap import (
    GAP_SIZE_ERROR_PERCENTILE,
    GAP_SIZE_WARN_PERCENTILE,
    check_gap_against_lengths,
    length_distance_thresholds,
    load_pickle,
)

data = load_pickle()
thr_warn = length_distance_thresholds(data, GAP_SIZE_WARN_PERCENTILE)
thr_error = length_distance_thresholds(data, GAP_SIZE_ERROR_PERCENTILE)
verdict = check_gap_against_lengths(41, 154, 174, thr_warn, thr_error)
```

Also exported: `load_pdb`, `find_gaps_in_chain`, `ca_distances`,
`aggregate_by_residue_range`, `lookup_gap_estimate`,
`pooled_terminal_budget`, `split_terminal_budget`, and
`smallest_feasible_res_min()` — the re-estimating scan behind the
length-range suggestions above. Pass `use_high=True` for the recommended
(upper-estimate) bound:

```python
from autocontigmap import load_pickle, smallest_feasible_res_min

data = load_pickle()
smallest_feasible_res_min(data, motif_residues=120, gap_angstroms=[40], res_max=400)
# 225
smallest_feasible_res_min(data, 120, [40], 400, use_high=True)
# 268
```
