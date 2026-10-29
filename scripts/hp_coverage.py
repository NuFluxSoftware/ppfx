#!/usr/bin/env python3
"""Python port of hp_coverage.C (PyROOT).

Compares two PPFX CreateIMapHists runs -- one over all interactions, one
with CUTNA49=1/CUTMIPP=1 (interactions not covered by thin-target/MIPP
hadron-production data) -- to quantify what fraction of interactions in a
given neutrino sample are actually backed by real hadron-production
measurements, for a given incident projectile species.

get_h2 (see ppfx_imap_common.py) is the only thing this macro used from the
external ~/root_utils/move_stats.C utility library. `center_titles()` is not
factored out as a generic utility here -- each call site only ever acts on
the one histogram just drawn in the current pad, so this script just calls
`h.GetXaxis().CenterTitle()`/`h.GetYaxis().CenterTitle()` directly.

See notes/ppfx-interaction-map.md and notes/ppfx-scripts-python-port.md.
"""
import argparse
import os

import ROOT

ROOT.gROOT.SetBatch(True)  # without this, batch PNG rendering subtly differs from `root -l -b -q`

from ppfx_imap_common import get_h2


def center_titles(h):
    h.GetXaxis().CenterTitle()
    h.GetYaxis().CenterTitle()


def style_matrix(h):
    h.SetTitle("; ;")
    h.LabelsDeflate("X")
    h.LabelsDeflate("Y")
    h.LabelsOption("a", "X")
    h.LabelsOption("a", "Y")
    h.SetMinimum(0.001)
    h.SetMarkerSize(1.5)
    h.SetMarkerColor(ROOT.kWhite)


def print_breakdown(h_all, h_nohp, row_label):
    print(f"{row_label:>13} | # interactions / event")
    print(f"{'particle':>13} |   nohp  |    all")
    print("-----------------------------------")
    for i in range(1, h_all.GetNbinsX() + 1):
        label = h_all.GetXaxis().GetBinLabel(i)
        print(f"{label:>13} |   {h_nohp.GetBinContent(i):3.3f} |   {h_all.GetBinContent(i):3.3f}")
    print("\n\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("file_all", help="CreateIMapHists output .root file, no coverage cuts")
    ap.add_argument("file_nohp", help="CreateIMapHists output .root file, CUTNA49=1 CUTMIPP=1")
    ap.add_argument("projectile", help='incident species as named in the imap file, e.g. "pi+"')
    ap.add_argument("--outdir", default=".", help="directory for output PNGs")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPaintTextFormat("3.3f")
    ROOT.gStyle.SetPadLeftMargin(0.2)
    ROOT.gStyle.SetPadBottomMargin(0.2)

    keep = []

    def save(c):
        c.SaveAs(os.path.join(args.outdir, f"hp_coverage_{c.GetName()}.png"))

    c_inter_all = ROOT.TCanvas("c_inter_all", "c_inter_all", 1000, 700)
    h_inter_all = get_h2("Projectile/h_in_vs_mat", args.file_all)
    nent = h_inter_all.GetEntries()
    style_matrix(h_inter_all)
    h_inter_all.Draw("col,text")
    h_inter_all.SetTitle("Interactions per event; target material; projectile")
    h_inter_all.GetYaxis().SetTitleOffset(1.8)
    h_inter_all.GetXaxis().SetTitleOffset(1.5)
    center_titles(h_inter_all)
    save(c_inter_all)

    c_inter_nohp = ROOT.TCanvas("c_inter_nohp", "c_inter_nohp", 1000, 700)
    h_inter_nohp = get_h2("Projectile/h_in_vs_mat", args.file_nohp)
    style_matrix(h_inter_nohp)
    h_inter_nohp.Draw("col,text")
    h_inter_nohp.SetTitle("Interactions per event unconstrained by HP data; target material; projectile")
    h_inter_nohp.GetYaxis().SetTitleOffset(1.8)
    h_inter_nohp.GetXaxis().SetTitleOffset(1.5)
    center_titles(h_inter_nohp)
    save(c_inter_nohp)

    c_proj = ROOT.TCanvas("c_proj", "c_proj", 1000, 700)
    h_proj_all = h_inter_all.ProjectionY("h_proj_all")
    h_proj_nohp = h_inter_nohp.ProjectionY("h_proj_nohp")
    h_proj_nohp.SetFillStyle(1001)
    h_proj_nohp.SetFillColor(ROOT.kGray)
    h_proj_all.LabelsOption("a", "X")
    h_proj_nohp.LabelsOption("a", "X")
    h_proj_all.Draw("hist")
    h_proj_all.SetTitle("; projectile ; interactions / event")
    h_proj_nohp.Draw("same")
    save(c_proj)

    n_all = h_proj_all.Integral()
    n_nohp = h_proj_nohp.Integral()

    print(f"Average # of interactions / event = {n_all}")
    print(f"              without HP coverage = {n_nohp}")
    print(f"                    HP coverage % = {(n_all - n_nohp) / n_all * 100}")
    print(f"Number of interactions = {nent}")
    print(f"Rough number of events = {nent / n_all}")
    print("\n\n")

    print_breakdown(h_proj_all, h_proj_nohp, "incident")

    c_projectile_all = ROOT.TCanvas("c_projectile_all", "c_projectile_all", 1000, 700)
    h_projectile_all = get_h2(f"Projectile/{args.projectile}/hmats{args.projectile}", args.file_all)
    style_matrix(h_projectile_all)
    h_projectile_all.Draw("col,text")
    h_projectile_all.SetTitle(f"{args.projectile} incident: Interactions per event; target material; produced particle")
    h_projectile_all.GetYaxis().SetTitleOffset(2.1)
    h_projectile_all.GetXaxis().SetTitleOffset(1.5)
    center_titles(h_projectile_all)
    save(c_projectile_all)

    c_projectile_nohp = ROOT.TCanvas("c_projectile_nohp", "c_projectile_nohp", 1000, 700)
    h_projectile_nohp = get_h2(f"Projectile/{args.projectile}/hmats{args.projectile}", args.file_nohp)
    style_matrix(h_projectile_nohp)
    h_projectile_nohp.Draw("col,text")
    h_projectile_nohp.SetTitle(
        f"{args.projectile} incident: Interactions per event unconstrained by HP data; "
        "target material; produced particle"
    )
    h_projectile_nohp.GetYaxis().SetTitleOffset(2.1)
    h_projectile_nohp.GetXaxis().SetTitleOffset(1.5)
    center_titles(h_projectile_nohp)
    save(c_projectile_nohp)

    c_prod = ROOT.TCanvas("c_prod", "c_prod", 1000, 700)
    h_prod_all = h_projectile_all.ProjectionY("h_prod_all")
    h_prod_nohp = h_projectile_nohp.ProjectionY("h_prod_nohp")
    h_prod_nohp.SetFillStyle(1001)
    h_prod_nohp.SetFillColor(ROOT.kGray)
    h_prod_all.LabelsOption("a", "X")
    h_prod_nohp.LabelsOption("a", "X")
    h_prod_all.Draw("hist")
    h_prod_all.SetTitle(f"{args.projectile} incident; produced particle ; interactions / event")
    h_prod_nohp.Draw("same")
    save(c_prod)

    print(f"{args.projectile:>13} incident")
    print_breakdown(h_prod_all, h_prod_nohp, "produced")

    keep += [
        h_inter_all, h_inter_nohp, h_proj_all, h_proj_nohp,
        h_projectile_all, h_projectile_nohp, h_prod_all, h_prod_nohp,
    ]


if __name__ == "__main__":
    main()
