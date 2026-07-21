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

Note: motif.pdb should only include the motif residues. Gaps are deteceted from that automatically. There is no way right now to input the whole protein and define the motif afterwards.

Output on stdout:

```
contig="5-15,60-100,5-15"
```

where the first and last entries are the N-/C-terminal tail budgets and the
middle entries are the estimated residue counts per internal gap.

Options:

- `--pickle-file NAME_OR_PATH` — which statistics checkpoint to use: one of
  the bundled variants `gyr` (default), `gyr_ss_2`, `standard`, `surface`,
  or a path to an external `.pkl` file.
- `--chain-order B,A` — for motifs whose segments are split across separate
  PDB chains (one chain per segment) instead of one chain with internal
  chain breaks. Gaps are then measured between the last residue of each
  chain and the first residue of the next.
- `--simple-terminals` — always emit `0-term_hi` at both ends instead of the
  default `res_min`/`res_max` conflict-clamped terminal budget.

## Python API

```python
from autocontigmap import contigmap

contig = contigmap("motif.pdb", res=[150, 200])
```

`contigmap()` prints its progress and returns the contig string (or `""` if
the requested length range can't accommodate the estimated gaps).
