# Containerizing ppfx

**Status: not yet built.** This document lays out a concrete path for
containerizing ppfx, based directly on a Docker/Singularity image that
*was* built, deployed to an HPC cluster (W&M's SciClone, SLURM +
Singularity), and verified working for the sibling package
[g4numi](https://github.com/NuSoftHEP/g4numi) on 2026-09-11 — see
`../../nusofthep/g4numi/docker/README.md` for that full, verified
writeup, and `~/obsidian/learn_root/notes/hpc-singularity-slurm.md` (a
separate personal notes vault, not part of either repo) for the session
history behind it. Nothing below has been built or run yet for ppfx
specifically; treat it as a starting plan for whoever picks this up
next, not as instructions known to work.

## Why this is a natural extension of the g4numi work

ppfx reads g4numi's own output format (`dk2nu` ntuples) as its primary
input — the two packages are meant to run back-to-back in a beamline
flux-prediction pipeline. They also share a dependency: both link
against `libdk2nuTree` via a `DK2NU` environment variable pointing at a
built dk2nu tree library (see `setup.sh`'s "SETUP ROOT, BOOST and DK2NU"
step, and g4numi's own `CMakeLists.txt` `WITH_DK2NU` option — same
library, same env var convention, in both packages). A container image
carrying Geant4 + dk2nu + g4numi (already built) is most of the way to
one that also carries ppfx — plausibly even the *same* image, or a
second image built `FROM` it, rather than a wholly separate build.

## The build-system complication, and the way around it

ppfx's primary/documented build path (`CMakeLists.txt`) depends on
Fermilab's `cetmodules`/`cetbuildtools` CMake extensions and a
`nufinder` package — both are UPS products from the Fermilab/NOvA
software ecosystem, not available as ordinary system or PyPI/conda
packages, and neither is present on this machine. Standing that up
inside a plain Docker image (outside UPS entirely) would be a real
undertaking on its own and hasn't been attempted.

**There is a second, much simpler path already in this repo**: a plain
top-level `Makefile`, explicitly fixed on 2026-08-17 (commit `f0a9775`,
"Fix plain Makefile build") to work correctly *outside* UPS — that
commit's own message describes fixing an unset `BOOST_INC` breaking the
include path, a link-order bug hiding `libppfx.so`'s undefined symbols,
and an unqualified `isinf` call that newer libstdc++/C++23 no longer
exposes globally. Reading that `Makefile` directly, its real
dependencies are just:

- `ROOT` (via `root-config --incdir`/`--libs`, plus `-lEG` for
  `TDatabasePDG`) — already the container's base image, same as g4numi.
- `Boost` headers only, no compiled Boost libraries (`BOOST_INC`,
  defaults to `/usr/include`) — a single `apt-get install libboost-dev`
  package away; already confirmed present on this dev machine
  (`libboost-dev` 1.90.0.1).
- `dk2nu` (`DK2NU` env var → `${DK2NU}/include` and
  `${DK2NU}/lib/libdk2nuTree.so`) — exactly the same dependency, built
  the same way, already solved in the g4numi image.

This is a much smaller and more container-friendly dependency set than
the `cetmodules` path, and reuses build logic already worked out for
g4numi almost directly.

## A concrete next-session plan

1. Start `FROM` the already-built `g4numi:11.4.2` image (or rebuild a
   shared base stage) — Geant4 isn't actually needed for ppfx itself,
   but ROOT and a working `DK2NU`-pointed dk2nu install already are, so
   reusing that layer avoids redoing work.
2. `apt-get install libboost-dev` (or use a base image that already has
   it).
3. `git archive HEAD` this repo the same way g4numi's Dockerfile does,
   to exclude any local build artifacts/output rather than a plain copy.
4. Run `make` (the plain top-level `Makefile`, **not** the
   `cetmodules`-based `CMakeLists.txt` path) with `DK2NU` pointed at the
   image's dk2nu install. This has not been tried inside a container —
   the 2026-08-17 Makefile fixes were verified on this dev machine
   directly, not in any container, so treat the first attempt as a real
   test, not a formality.
5. Verify with one of the small executables under `src/` (e.g.
   `PrintMaterials.C`, which has minimal input requirements) before
   attempting a full reweighting run (`doReweight_dk2nu_originalpip.C`),
   which needs a real g4numi `dk2nu` output file as input — conveniently,
   exactly the kind of file the g4numi image can now produce, making an
   end-to-end g4numi→ppfx containerized pipeline test a natural
   follow-up once both halves build.
6. Deployment to an HPC cluster (convert to `.sif`, stage via `scp`,
   `--bind`-mount for I/O, submit any heavy `singularity build`/`pull`
   step as a batch job rather than running it on a login node) follows
   exactly the pattern already verified for g4numi — see that package's
   `docker/README.md`, nothing ppfx-specific expected to differ there.
