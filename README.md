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

A motif with no gaps at all (a single chain via `--chain-order`) doesn't
split the budget: each terminal gets `0` to the whole upper budget,
`res_max - motif`, and the `length` field pins the total.

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
  Error: res_min (200) is below the shortest buildable design (17 aa motif + 184-318 aa of gaps) -- nearest buildable length: 201 (or 199 below); or lower the lower quantile: --quantiles 0.25 0.6 builds 200.
  ```

  A `res_max` below the same floor reports the same error — `res_min <=
  res_max`, so a `res_max` under the floor puts `res_min` under it too.

  The message offers two ways out:

  - **Nearest buildable length** above and below. For a range the whole
    window moves and keeps its width (`nearest buildable range: 192-202 (or 189-199 below)`).
    Every candidate **re-estimates the gaps for its own window**: sequence
    separation changes with chain length, so the `motif + gap_lo_sum` floor
    computed for the window that just failed is not itself a usable length.
    That also makes the answer non-monotonic — a shorter design can build
    when a longer one does not, because the gap estimates shrink too.
  - **A lower lower quantile** (`--quantiles`) that makes the requested
    length build with the upper quantile kept — the highest stored one that
    works. If none down to 0.25 does, the message says so; if you are
    already at 0.25 there is nothing to suggest.

- **Warning** (not fatal) if the lower terminal budget is negative but
  every requested length is still reachable with the gaps near their lower
  estimates — the budget is pinned to 0 and the contig is still printed:

  ```
  Warning: the lower terminal budget is negative and was pinned to 0; res_min (250) is still reachable, but only when the gaps sample near their lower estimates (109 aa total). Nearest length without pinning: <above> (or <below> below).
  ```

  The suggestion uses the same nearest-length scan as the error above, but
  asks for the *high* end of every gap estimate to fit, which is exactly the
  threshold between a pinned and an unpinned lower budget. The tip is
  dropped if no length in the checkpoint range clears it.

Options:

- `--pickle-file NAME_OR_PATH` — which statistics checkpoint to use. All
  bundled variants come from non-redundant PISCES chains (40% sequence
  identity, X-ray, <= 3.0 A):
  - `pisces_pc40_s5_nohelix` (default) — residue pairs at least 5 apart in
    sequence (gaps between motif segments rather than residues within one
    segment), minus pairs whose whole stretch lies in one continuous
    alpha-helix (flagged from Calpha distances: d(i,i+3) 4.6-5.8 A and
    d(i,i+4) 5.6-6.8 A). Without this, intra-helix i->i+5..i+8 pairs pile up
    at 8-11 A and pull the gap estimates there down by up to ~9 residues.
    Quantiles are pooled over all chains.
  - `pisces_pc40_s5` — pairs at least 5 apart in sequence, helices kept.
  - `pisces_pc40` — all pair separations.
  - or a path to an external `.pkl` file.
- `--quantiles LOW HIGH` — which pair of checkpoint quantiles bounds each
  gap estimate (default `0.4 0.6`). `--quantiles 0.25 0.75` gives every gap a
  wider range, which makes a fixed design length easier to fit when many gaps
  are summed. Available: 0.25, 0.4, 0.5, 0.6, 0.75. In the library:
  `build_contig(..., quantiles=quantile_keys(0.25, 0.75))`.
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
covered `[40, 499]` range and raises `ValueError` if the requested range
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
`pooled_terminal_budget`, `split_terminal_budget`,
`ungapped_terminal_budget`, and
and the two helpers behind the length-error suggestions above:
`nearest_feasible_lengths()` returns `(below, above)`, the nearest
buildable `res_min` on each side with the window width kept (`use_high=True`
for "no pinning"), and `lower_quantile_fix()` returns `(lower_key, tried)`:

```python
from autocontigmap import load_pickle, nearest_feasible_lengths, lower_quantile_fix

data = load_pickle()
nearest_feasible_lengths(data, motif_residues=120, gap_angstroms=[40], res_min=181, res_max=181)
lower_quantile_fix(data, 120, [40], 181, 181)   # ("Q0.25", True) or (None, ...)
```
