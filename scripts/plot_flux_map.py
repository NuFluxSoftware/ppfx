#!/usr/bin/env python3
"""Python port of plot_flux_map.C (PyROOT).

Plots PPFX CreateIMapHists interaction-map output for one incident particle
species: the overall material x projectile interaction matrix, that
species' material/volume breakdowns, and its incident kinetic-energy
spectrum (linear + log).

get_h1/get_h2 (see ppfx_imap_common.py) are the only two things this macro
used from the external ~/root_utils/move_stats.C utility library.

See notes/ppfx-interaction-map.md and notes/ppfx-scripts-python-port.md.
"""
import argparse
import os

import ROOT

ROOT.gROOT.SetBatch(True)  # without this, batch PNG rendering subtly differs from `root -l -b -q`

from ppfx_imap_common import get_h1, get_h2


def slugify(particle):
    return particle.replace("+", "p").replace("-", "m").replace(" ", "")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("particle", help='incident species as named in the imap file, e.g. "pi+"')
    ap.add_argument("file", help="CreateIMapHists output .root file")
    ap.add_argument("--outdir", default=".", help="directory for output PNGs")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat("3.3f")

    slug = slugify(args.particle)
    keep = []

    def save(c):
        c.SaveAs(os.path.join(args.outdir, f"flux_map_{c.GetName()}_{slug}.png"))

    cinter = ROOT.TCanvas("cinter", "cinter", 1000, 700)
    hinter = get_h2("Projectile/h_in_vs_mat", args.file)
    hinter.SetTitle("; ;")
    hinter.LabelsDeflate("X")
    hinter.LabelsDeflate("Y")
    hinter.LabelsOption("a", "X")
    hinter.LabelsOption("a", "Y")
    hinter.SetMinimum(0.001)
    hinter.SetMarkerSize(1.5)
    hinter.SetMarkerColor(ROOT.kWhite)
    hinter.Draw("col,text")
    save(cinter)
    keep.append((cinter, hinter))

    cmat = ROOT.TCanvas("cmat", "cmat", 1000, 700)
    ROOT.gPad.SetLeftMargin(0.15)
    h_mats = get_h2(f"Projectile/{args.particle}/hmats{args.particle}", args.file)
    h_mats.LabelsDeflate("X")
    h_mats.LabelsDeflate("Y")
    h_mats.LabelsOption("a", "X")
    h_mats.LabelsOption("a", "Y")
    h_mats.SetMinimum(0.001)
    h_mats.SetMarkerSize(1.5)
    h_mats.SetMarkerColor(ROOT.kWhite)
    h_mats.Draw("col,text")
    save(cmat)
    keep.append((cmat, h_mats))

    cvol = ROOT.TCanvas("cvol", "cvol", 1000, 700)
    ROOT.gPad.SetLeftMargin(0.15)
    h_vols = get_h2(f"Projectile/{args.particle}/hvols{args.particle}", args.file)
    h_vols.LabelsDeflate("X")
    h_vols.LabelsDeflate("Y")
    h_vols.LabelsOption("a", "Y")
    h_vols.SetMinimum(0.001)
    h_vols.SetMarkerSize(1.5)
    h_vols.SetMarkerColor(ROOT.kWhite)
    h_vols.Draw("col,text")
    save(cvol)
    keep.append((cvol, h_vols))

    cine = ROOT.TCanvas("cine", "cine", 1400, 700)
    h_ine = get_h1(f"Projectile/{args.particle}/henergy{args.particle}", args.file)
    h_ine.Rebin(10)
    cine.Divide(2, 1)
    cine.cd(1)
    h_ine.SetTitle(f"Incident {args.particle} energy; kinetic energy (GeV); ")
    h_ine.SetFillColor(ROOT.kOrange - 3)
    h_ine.GetXaxis().SetRangeUser(0, 125)
    h_ine.Draw("hist")
    h_ine.GetXaxis().SetNoExponent()
    cine.cd(2)
    ROOT.gPad.SetLogy()
    hc_ine = h_ine.DrawClone()
    hc_ine.SetTitle("Same, in logscale")
    cine.cd(0)
    save(cine)
    keep.append((cine, h_ine, hc_ine))

    print("Done.")
    return keep  # kept alive for interactive use; CLI run exits right after saving


if __name__ == "__main__":
    main()
