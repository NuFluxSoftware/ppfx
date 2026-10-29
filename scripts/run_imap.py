#!/usr/bin/env python3
"""Thin CLI wrapper around PPFX's CreateIMapHists binary.

Replaces the current ad hoc invocation pattern (the relevant one line out of
$PPFX_DIR/scripts/imap_job.sh, which otherwise targets Fermilab's FIFE grid
and isn't usable on this machine -- see notes/ppfx-interaction-map.md): sets
up PPFX_DIR/DK2NU/MODE/ROOTSYS and the right LD_LIBRARY_PATH automatically,
runs the program, and prints a one-line summary (total weight, warning
count) pulled from its log instead of leaving that to manual grepping.

Works against either PPFX checkout via --ppfx-dir; defaults to
~/ppfx_dev/nufluxsoftware/ppfx, the current development target (see
notes/ppfx-nufluxsoftware-checkout.md) -- pass --ppfx-dir ~/minerva/ppfx to
run the same job against the older, frozen fork for comparison.

Example (reproduces notes/ppfx-interaction-map.md's verification run):
  run_imap.py --elow 6 --ehigh 8 --nu 14 --outfile all_6_8.root \\
      --inputs "inputs/*.root"
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
    ap.add_argument("--inputs", required=True, help='input flux file glob, e.g. "inputs/*.root" (passed through to TChain::Add, not expanded by this wrapper)')
    ap.add_argument("--elow", type=float, required=True, help="lower edge of the neutrino energy window (GeV)")
    ap.add_argument("--ehigh", type=float, required=True, help="upper edge of the neutrino energy window (GeV)")
    ap.add_argument("--nu", type=int, default=14, help="neutrino PDG code (default: 14, numu)")
    ap.add_argument("--cut-na49", type=int, default=0, choices=[0, 1], help="skip interactions covered by NA49 thin-target data")
    ap.add_argument("--cut-mipp", type=int, default=0, choices=[0, 1], help="skip interactions covered by MIPP data")
    ap.add_argument("--outfile", required=True, help="output .root filename")
    ap.add_argument("--logfile", default=None, help="log filename (default: outfile with .root -> .log)")
    ap.add_argument("--cwd", default=".", help="directory to run in (inputs glob is resolved relative to this)")
    args = ap.parse_args()

    binary = os.path.join(args.ppfx_dir, "bin", "CreateIMapHists")
    if not os.path.isdir(args.ppfx_dir):
        sys.exit(f"--ppfx-dir does not exist: {args.ppfx_dir}")
    if not os.path.isfile(binary):
        sys.exit(f"CreateIMapHists not found (build it first?): {binary}")
    if not os.path.isdir(args.dk2nu):
        sys.exit(f"--dk2nu does not exist: {args.dk2nu}")

    logfile = args.logfile or re.sub(r"\.root$", "", args.outfile) + ".log"

    env = dict(os.environ)
    env["ROOTSYS"] = args.rootsys
    env["PATH"] = f"{args.rootsys}/bin:" + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = (
        f"{args.rootsys}/lib:{args.dk2nu}/lib:{args.ppfx_dir}/lib:" + env.get("LD_LIBRARY_PATH", "")
    )
    env["DK2NU"] = args.dk2nu
    env["PPFX_DIR"] = args.ppfx_dir
    env["MODE"] = args.mode

    argv = [
        binary, args.outfile, args.inputs,
        str(args.elow), str(args.ehigh), str(args.nu),
        str(args.cut_na49), str(args.cut_mipp),
    ]
    print(f"+ {' '.join(argv)}  (PPFX_DIR={args.ppfx_dir})")

    with open(os.path.join(args.cwd, logfile), "w") as lf:
        result = subprocess.run(argv, cwd=args.cwd, env=env, stdout=lf, stderr=subprocess.STDOUT)

    # errors="replace": ROOT's own console output occasionally contains stray
    # non-UTF-8 bytes (e.g. progress-bar control characters); this is just a
    # log we're grepping for summary lines, not something that needs exact
    # byte fidelity.
    log_text = open(os.path.join(args.cwd, logfile), errors="replace").read()
    total_weight = re.search(r"Total weight:\s*([\d.eE+-]+)", log_text)
    n_skipped = len(re.findall(r"Skipping unknown volume", log_text))

    print(f"exit code: {result.returncode}  (log: {logfile})")
    if total_weight:
        print(f"Total weight: {total_weight.group(1)}")
    print(f'"Skipping unknown volume" warnings: {n_skipped}')

    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
