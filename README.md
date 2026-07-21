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

Output on stdout:

```
contig="5-15,60-100,5-15"
```

where the first and last entries are the N-/C-terminal tail budgets and the
middle entries are the estimated residue counts per internal gap.

Pass `--rfd1` to instead print a full, old-style RFdiffusion(1) contig
string, with chain letters and fixed residue ranges included:

```
contigmap.contigs=[5-15/B165-178/60-100/A0-20/5-15]
```

Options:

- `--pickle-file NAME_OR_PATH` — which statistics checkpoint to use: one of
  the bundled variants `gyr` (default), `gyr_ss_2`, `standard`, `surface`,
  or a path to an external `.pkl` file.
- `--rfd1` — output the old-style `contigmap.contigs=[...]` contig instead
  of the plain numeric-ranges-only one. See "Chain selection" below for how
  multi-chain PDBs are handled in this mode.
- `--chain-order B,A` — for motifs whose segments are split across separate
  PDB chains (one chain per segment) instead of one chain with internal
  chain breaks. Gaps are then measured between the last residue of each
  chain and the first residue of the next.
- `--simple-terminals` — always emit `0-term_hi` at both ends instead of the
  default `res_min`/`res_max` conflict-clamped terminal budget.

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
a gap are carried through unchanged as fixed spans. In `--rfd1` output, all
of these segments are joined with a hard chain break (`/0 `) in the PDB's
chain order — each ends up as its own separate output chain.

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
