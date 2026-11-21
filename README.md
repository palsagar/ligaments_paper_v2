# Drop-size statistics of randomly corrugated ligaments

This repository holds the code, data pipeline and LaTeX source for a follow-up
to the paper *Statistics of drops generated from ensembles of randomly
corrugated ligaments*, Phys. Fluids 36, 112116 (2024),
[doi:10.1063/5.0221732](https://doi.org/10.1063/5.0221732).

The aim is a theory that predicts the drop-size distribution of a liquid
ligament from its initial surface noise, its noise amplitude and its Ohnesorge
number, with no parameter fitted to the simulated ensembles, and a simulation
campaign that tests that theory on parameter sets it has not seen.

**Status: set-up.** The solver build and one test case work. There are no
results yet.

## Layout

| Path | Content |
| --- | --- |
| `sim/` | Basilisk C cases. `ligament.c` is a prototype: one axisymmetric, periodic, randomly corrugated ligament. |
| `theory/` | Python models. `linear_check.py` evolves random surfaces with linear Rayleigh–Plateau theory. |
| `paper/` | LaTeX source of the paper. |
| `third_party/` | The exact Basilisk source tarball used, with its SHA-256 checksum. |
| `Makefile` | Entry point for every step. |

## Requirements

- A C99 compiler and GNU make. On Linux also `gawk`.
- [`uv`](https://docs.astral.sh/uv/) for the Python code.
- `latexmk` and a TeX distribution for the paper.

## Use

```sh
make basilisk       # build the Basilisk compiler (qcc) from the archived tarball
make smoke          # run one coarse ligament and check it (about 15 s)
make linear-check   # linear-theory model on 20,000 random surfaces (about 3 min)
make paper          # build paper/build/main.pdf
```

## The solver is pinned

The simulations use [Basilisk](https://basilisk.fr/). Its download URL always
serves the latest version, so the tarball used here is archived in
`third_party/` (files dated 2026-10-25, release tag "20-10-26") and the build
verifies its checksum.

## Licences

- Code: GNU General Public License v3 (`LICENSE`). Basilisk is GPL, and the
  simulation cases include its headers.
- Paper text and data: Creative Commons Attribution 4.0
  (`LICENSE-CC-BY-4.0`).
- `third_party/basilisk-2026-10-25.tar.gz` is redistributed under its own
  licence, the GPL (see `basilisk/src/COPYING` inside the tarball).
