# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

A research project that produces one journal paper: a follow-up to *Statistics of drops generated from ensembles of randomly corrugated ligaments* (Phys. Fluids 36, 112116, 2024). The goal is a theory that predicts the drop-size distribution of a liquid ligament from its initial surface noise, noise amplitude η and Ohnesorge number Oh, with no parameter fitted to simulated ensembles, tested on parameter sets simulated only after the predictions are committed. The repository is public and must rebuild the paper from raw data.

The full plan (goal, parts, pass conditions, tested assumptions, what was left out) is in `docs/deeper/2026-11-20-ligaments-paper-v2.md`. Read it before designing anything. That folder is untracked on purpose, see "Files that must not be committed".

## Commands

```sh
make basilisk       # build qcc from third_party/basilisk-2026-10-25.tar.gz into build/basilisk/src (verifies the SHA-256)
make smoke          # compile sim/ligament.c and run one ligament at level 9, seed 1 (about 15 s)
make linear-check   # uv run python theory/linear_check.py 20000 400 (about 3 min)
make paper          # latexmk, output in paper/build/
make clean
```

Run one ligament by hand after `make smoke` has built the binary:

```sh
cd build/smoke && ../ligament MAXLEVEL SEED [POISSON_TOLERANCE]
```

Compile any other Basilisk case the same way the Makefile does (`BASILISK` must point at `build/basilisk/src`):

```sh
BASILISK=$PWD/build/basilisk/src build/basilisk/src/qcc -O2 -Wall case.c -o case -lm
```

There is no unit-test suite and no linter yet. `make smoke` is the only automated check: it passes when the run reaches the end, the thread has broken, and the liquid volume is conserved to 0.1%. CI (`.github/workflows/ci.yml`) runs `make basilisk`, `make smoke` and a LaTeX build on Ubuntu for every push.

Python runs through `uv` only (`uv run ...`, `uv add ...`). LaTeX is not installed on the development Mac; the paper is built in CI.

## Architecture

Three kinds of code, joined by the Makefile:

- **`sim/`: Basilisk C cases.** Basilisk C is a C dialect; `qcc` translates it to C99 and compiles it. Cases are axisymmetric (`axi.h`), two-phase with surface tension, on an adaptive quadtree, periodic along the axis. `ligament.c` is a prototype used for timing and the smoke test, not the production case.
- **`theory/`: Python models** of the same physics at lower cost. `linear_check.py` generates random surfaces, evolves each Fourier mode with the linear Rayleigh–Plateau dispersion relation, and pinches where the linear radius reaches zero.
- **`paper/`: LaTeX.** Every number in the text is to be written by the analysis as a LaTeX macro, with one script per figure.

The solver is pinned by archiving its tarball in `third_party/`, because basilisk.fr serves only the latest version. Do not replace the tarball or build against a system-wide Basilisk.

## Conventions in the simulation code

- Units: mean ligament radius R = 1, liquid density 1, surface tension 1, so the capillary time is 1. `W = 2R` is the diameter and drop sizes are reported as `D/W`. Oh = liquid viscosity in these units.
- With `axi.h`, x is the axis, y is the radius, and `dv()` is the volume per radian. Multiply by 2π for a true volume.
- Mesh adaptation must act on a copy of the volume fraction with bilinear prolongation (`-DFULL_INTERFACE=1` in `ligament.c`). Adapting on `f` itself lets interface cells fall below the maximum level and the liquid volume drift by 0.5%.
- One serial process per ligament. OpenMP is about 7 times slower for this case, and Basilisk's GPU grids do not support quadtrees or surface tension.
- Refinement level N means a cell size of 100 R / 2^N for the aspect-ratio-50 ligament: level 10 is about 10 cells per radius. Measured cost to T = 14 on an M2 core: 77 s at level 10, about 400 s at level 11.
- Initial surface profiles for production runs come from Python as files of Fourier coefficients, never from C's `rand()`, so a seed gives the same ligament on every machine.
- A drop is a liquid body with D/W ≥ 0.25. Smaller bodies are debris and are logged separately.

## Rules of the study

These are scientific constraints, not style. Breaking one invalidates the paper's claim.

- **No fitting to ensembles.** The theory may use constants and tables only from deterministic single-object simulations (one wavelength of a perturbed cylinder, one retracting filament, single-fragment fate tables).
- **Open and closed parameter sets.** Only the baseline (η = 0.08, Oh = 0.01, 25 modes of equal power, aspect ratio 50) may be explored freely. Every other set is closed: its predictions and pass conditions are committed and pushed before it is simulated, and it is run once. Half of the baseline seeds are also held back.
- **The ligament is the statistical unit.** Drops from one ligament share a fixed volume. Resample whole ligaments for every error bar; never resample individual drops, and never fit histogram bins.
- **Record events, not only snapshots.** Production runs log every pinch-off and merger with time, position, volumes and velocities.

## Files that must not be committed

`.gitignore` excludes them; keep it that way unless the author says otherwise.

- `ligaments_paper_v1.pdf` (publisher's PDF) and `JFM_Perspective_first_revision-9.pdf` (an unpublished draft by other authors; nothing from it may be quoted or sent to outside services).
- `v1_review.md`.
- `docs/deeper/` (the plan and its evidence), private until the author decides to publish it.
