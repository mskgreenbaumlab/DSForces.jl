# DSForces

`dsforces` is a Python package to compute double-stranded (DS) forces in
nucleic-acid sequences.

This is a dependency-free Python port of the original Julia package. It keeps
Julia-compatible exported names such as `ComputeDSForce`, while also exposing
Python-style aliases such as `compute_ds_force`.

## Installation

From this repository:

```bash
python -m pip install -e .
```

## Usage

```python
from dsforces import ComputeDSForce

seq = "ACAACGTAACGGTCGAGTCG"

ds_force = ComputeDSForce(seq)
ds_force, seq_a, seq_b = ComputeDSForce(seq, return_lcs_positions=True)
```

Returned coordinates are 1-based inclusive `(start, end)` tuples, matching the
genomic coordinate convention used by the Julia package.

To constrain one of the complementary subsequences to a region, pass either a
1-based inclusive pair or a normal Python `range`:

```python
ComputeDSForce(seq, return_lcs_positions=True, seq_a_range=(1, 10))
ComputeDSForce(seq, return_lcs_positions=True, seq_a_range=range(1, 11))
```

Sliding-window computation uses stride 1:

```python
ComputeDSForce(seq * 2, return_lcs_positions=True, sliding_window_length=len(seq))
```

Ambiguous bases such as `N` return `None`, mirroring Julia's `missing`.

## Command line

The genome scanning script is available as:

```bash
scan-genome-dsforce genome.fa chr1 199000 201999 --window_length 3000 --outfile chr1.csv
```

## Tests

```bash
python -m unittest discover -s tests
```
