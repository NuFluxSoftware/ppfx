#!/usr/bin/env python3
"""Thin CLI wrapper around PPFX's reweighting ("rwgt") binary.

Auto-detects which binary this checkout has -- `doReweight_dk2nu` (3 args:
input, output, options -- ~/minerva/ppfx) or `doReweight_dk2nu_originalpip`
(4 args with a precalculated detector index, or 6 with an explicit x/y/z
position -- ~/ppfx_dev/nufluxsoftware/ppfx) -- and builds the right argv, so
one tool works against either checkout. Sets up PPFX_DIR/DK2NU/MODE/ROOTSYS
and the right LD_LIBRARY_PATH automatically.

Defaults to ~/ppfx_dev/nufluxsoftware/ppfx, the current development target
(see notes/ppfx-nufluxsoftware-checkout.md) -- pass --ppfx-dir ~/minerva/ppfx
to run the same job against the older, frozen fork for comparison (this is
exactly how that note's rwgt cross-check was done).

Example (reproduces notes/ppfx-weight-calculations.md's verification run):
  run_rwgt.py --input g4numiv6_dk2nu_minervame_me000z200i_0000_0006.root \\
      --output ppfx_out_0000.root
"""
import argparse
import os
import re
import subprocess
import sys


def default_ppfx_dir():
    return os.environ.get("PPFX_DIR", os.path.expanduser("~/ppfx_dev/nufluxsoftware/ppfx"))


def default_dk2nu():
    return os.environ.get("DK2NU", os.path.expanduser("~/ppfx_dev/nusofthep/dk2nu-build"))


def default_rootsys():
    return os.environ.get("ROOTSYS", os.path.expanduser("~/root"))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppfx-dir", default=default_ppfx_dir(), help="PPFX source tree (default: $PPFX_DIR or ~/ppfx_dev/nufluxsoftware/ppfx)")
    ap.add_argument("--dk2nu", default=default_dk2nu(), help="dk2nu install prefix (default: $DK2NU or ~/ppfx_dev/nusofthep/dk2nu-build)")
    ap.add_argument("--rootsys", default=default_rootsys(), help="ROOT install prefix (default: $ROOTSYS or ~/root)")
    ap.add_argument("--mode", default="NUMI", help="beamline mode (default: NUMI)")
    ap.add_argument("--input", required=True, help="input g4numi dk2nu/dkmeta flux file")
    ap.add_argument("--output", required=True, help="output .root filename")
    ap.add_argument("--options", default=None, help="PPFX options XML (default: $PPFX_DIR/scripts/inputs_default.xml)")
    ap.add_argument("--idet", default="1", help='precalculated detector index (doReweight_dk2nu_originalpip only; default "1", MINERvA)')
    ap.add_argument("--pos", nargs=3, metavar=("XPOS", "YPOS", "ZPOS"), default=None, help="explicit x/y/z position in cm (doReweight_dk2nu_originalpip only; overrides --idet)")
    ap.add_argument("--logfile", default=None, help="log filename (default: output with .root -> .log)")
    ap.add_argument("--cwd", default=".", help="directory to run in")
    args = ap.parse_args()

    options = args.options or os.path.join(args.ppfx_dir, "scripts", "inputs_default.xml")
    if not os.path.isdir(args.ppfx_dir):
        sys.exit(f"--ppfx-dir does not exist: {args.ppfx_dir}")
    if not os.path.isdir(args.dk2nu):
        sys.exit(f"--dk2nu does not exist: {args.dk2nu}")
    if not os.path.isfile(options):
        sys.exit(f"options XML not found: {options}")

    bindir = os.path.join(args.ppfx_dir, "bin")
    originalpip = os.path.join(bindir, "doReweight_dk2nu_originalpip")
    plain = os.path.join(bindir, "doReweight_dk2nu")
    if os.path.isfile(originalpip):
        binary = originalpip
        if args.pos:
            argv = [binary, args.input, args.output, options] + list(args.pos)
        else:
            argv = [binary, args.input, args.output, options, args.idet]
    elif os.path.isfile(plain):
        binary = plain
        argv = [binary, args.input, args.output, options]
    else:
        sys.exit(f"no rwgt binary found in {bindir} (build it first?)")

    logfile = args.logfile or re.sub(r"\.root$", "", args.output) + ".log"

    env = dict(os.environ)
    env["ROOTSYS"] = args.rootsys
    env["PATH"] = f"{args.rootsys}/bin:" + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = (
        f"{args.rootsys}/lib:{args.dk2nu}/lib:{args.ppfx_dir}/lib:" + env.get("LD_LIBRARY_PATH", "")
    )
    env["DK2NU"] = args.dk2nu
    env["PPFX_DIR"] = args.ppfx_dir
    env["MODE"] = args.mode

    print(f"+ {' '.join(argv)}  (PPFX_DIR={args.ppfx_dir})")

    with open(os.path.join(args.cwd, logfile), "w") as lf:
        result = subprocess.run(argv, cwd=args.cwd, env=env, stdout=lf, stderr=subprocess.STDOUT)

    # errors="replace": ROOT's own console output occasionally contains stray
    # non-UTF-8 bytes (e.g. progress-bar control characters); this is just a
    # log we're grepping for summary lines, not something that needs exact
    # byte fidelity.
    log_text = open(os.path.join(args.cwd, logfile), errors="replace").read()
    n_entries = re.search(r"N of entries:\s*(\d+)", log_text)
    finished = "End of run()" in log_text

    print(f"exit code: {result.returncode}  (log: {logfile})")
    if n_entries:
        print(f"N of entries: {n_entries.group(1)}")
    print(f"completed ('End of run()' seen in log): {finished}")

    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
