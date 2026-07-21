# AutoContigmap

Estimate RFdiffusion contig gap sizes for a motif PDB, from real PDB
residue-count-vs-Cα-distance statistics, instead of guessing a gap/length
range by hand.

Given a motif PDB with one or more chain breaks (or several single-segment
chains meant to be scaffolded together) and a target total protein-length
range, it looks up each gap's Cα-Cα distance against a bundled statistics
checkpoint and reports a plausible residue-count range per gap, plus N-/C-
terminal tail budgets — ready to paste into an RFdiffusion contig string.

## Install

```bash
pip install "git+https://github.com/lorenzkleiter/AutoContigmap.git"
```

## CLI usage

```bash
autocontigmap motif.pdb <res_min> <res_max>
```

Note: motif.pdb should only include the motif residues. Gaps are detected from that automatically. There is no way to input the whole protein and define the motif afterwards.

Output on stdout, with fixed motif/target chain spans included alongside
the estimated gap-fill ranges:

```
"contig": "8-15,A18-25,16-30,A47-54,16-30,A92-99"
"length": "150-200"
```

`length` is just `res_min-res_max` echoed back, for pinning the overall
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

- **Warning** (not fatal) if `res_min` clears that bar but `res_max` is
  still too tight for even the smallest per-gap estimate — the terminal
  budget clamps to 0 and the contig is still printed, but it likely won't
  satisfy `res_max`:

  ```
  Warning: even the smallest gap estimates (24 aa total) push the minimum feasible length to 56 aa, above res_max (40) -- a bigger design length range is probably necessary (res_max >= 56).
  ```

Options:

- `--pickle-file NAME_OR_PATH` — which statistics checkpoint to use: one of
  the bundled variants `gyr` (default), `gyr_ss_2`, `standard`, `surface`,
  or a path to an external `.pkl` file.
- `--rfd1` — output the old-style `contigmap.contigs=[...]`/`contigmap.length=...`
  form instead of the standard `"contig": ...`/`"length": ...` one.
- `--chain-order B,A` — for motifs whose segments are split across separate
  PDB chains (one chain per segment) instead of one chain with internal
  chain breaks. Gaps are then measured between the last residue of each
  chain and the first residue of the next.
- `--strict-terminals` — use the `res_min`/`res_max` conflict-clamped
  terminal budget instead of the default simple `0-term_hi` one.

## Chain selection

Any chain that contains at least one internal gap (chain break) is treated
as a motif chain and gap-filled. If no chain has a gap and `--chain-order`
isn't given, this is an error — there's nothing to scaffold:

```
Error: Not motif chain identified, motif chains have to have gaps. If the
motif is defined over multiple chains use --chain-order
```

**Default mode (no `--chain-order`)**: every chain with an internal gap is
gap-filled independently, each with its own terminal budget. Chains without
a gap are carried through unchanged as fixed spans. All of these segments
are joined with a hard chain break (`/0 `) in the PDB's chain order — each
ends up as its own separate output chain.

**`--chain-order` mode**: the listed chains are instead scaffolded into a
single continuous designed chain (no break between them, one shared
terminal budget), with gaps measured between consecutive chains. Any other
chains present in the PDB are still carried through unchanged as fixed
spans, each behind its own chain break.

## Python API

```python
from autocontigmap import contigmap

contig = contigmap("motif.pdb", res=[150, 200])
```

`contigmap()` prints its progress and returns the contig string (or `""` if
the requested length range can't accommodate the estimated gaps).

For callers that already have their own Cα-Cα distance measurement (not a
Bio.PDB structure) and just want the lookup, use the lower-level functions
directly — this is what the CLI itself is built from:

```python
from autocontigmap import estimate_gap_fill, gap_size_percentile_threshold, load_pickle

data = load_pickle()  # load once, reuse across calls
aa_low, aa_high = estimate_gap_fill(gap_ang=41, res_min=154, res_max=174, gap_size_data=data)
p95 = gap_size_percentile_threshold(res_min=154, res_max=174, gap_size_data=data)
```

Both clamp `res_min`/`res_max` to the checkpoint's covered `[10, 499]`
range and raise `ValueError` if the requested range doesn't overlap it at
all.
