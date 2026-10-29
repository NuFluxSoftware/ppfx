#!/usr/bin/env python3
"""Python port of plot_ppfx_ratio.C (PyROOT).

Weighted/unweighted (CV-reweighted / nominal) flux ratio vs. energy, with a
many-universes 1-sigma band computed directly on the ratio (each universe's
own ratio-to-nominal, not propagated from the absolute band).

See notes/ppfx-weight-calculations.md and notes/ppfx-scripts-python-port.md.
"""
import argparse
import os

import ROOT

ROOT.gROOT.SetBatch(True)  # without this, batch PNG rendering subtly differs from `root -l -b -q`

from ppfx_common import COLORS, FLAVOR_LABELS, FLAVORS, universe_band


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="ppfx_out_0000.root", help="PPFX doReweight_dk2nu output file")
    ap.add_argument("--outdir", default=".", help="directory for output PNG")
    ap.add_argument("--nuniv", type=int, default=100, help="number of PPFX universes in the file")
    ap.add_argument("--emax", type=float, default=60.0, help="upper edge of the energy frame (GeV)")
    ap.add_argument(
        "--label",
        default="g4numiv6 2018 production, me000z200i (ME beam, +200kA horn), file 0000",
        help="bottom annotation text",
    )
    ap.add_argument("--output", default="ppfx_ratio_band.png", help="output PNG filename")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(0.13)
    ROOT.gStyle.SetPadBottomMargin(0.12)

    f = ROOT.TFile(args.input)

    c1 = ROOT.TCanvas("c1", "c1", 1100, 850)
    c1.Divide(2, 2)

    # keep every ROOT object alive for the lifetime of the canvas (PyROOT has
    # no owning C++ stack frame to do this for us the way the macro's
    # function-local variables did)
    keep = []

    for i, flav in enumerate(FLAVORS):
        hnom = f.Get(f"nom/hnom_{flav}")
        hcv = f.Get(f"nom/hcv_{flav}")
        hnom.SetDirectory(0)
        hcv.SetDirectory(0)

        huniv = [f.Get(f"{flav}_total/htotal_{flav}_{j}") for j in range(args.nuniv)]

        nbins = hnom.GetNbinsX()

        gcv = ROOT.TGraph()
        npcv = 0
        gband = ROOT.TGraphAsymmErrors()
        npband = 0

        for b in range(1, nbins + 1):
            nomv = hnom.GetBinContent(b)
            if nomv <= 0:
                continue  # undefined ratio, skip (matches earlier plots' single-file-stats caveat)

            x = hnom.GetBinCenter(b)
            xw = hnom.GetBinWidth(b) / 2.0

            cvv = hcv.GetBinContent(b)
            gcv.SetPoint(npcv, x, cvv / nomv)
            npcv += 1

            ratios = [h.GetBinContent(b) / nomv for h in huniv]
            mean = sum(ratios) / len(ratios)
            mean2 = sum(r * r for r in ratios) / len(ratios)
            rms = max(0.0, mean2 - mean * mean) ** 0.5

            gband.SetPoint(npband, x, mean)
            gband.SetPointError(npband, xw, xw, rms, rms)
            npband += 1

        gband.SetFillColorAlpha(COLORS[i], 0.35)
        gband.SetLineColor(0)

        gcv.SetLineColor(COLORS[i])
        gcv.SetLineWidth(2)
        gcv.SetMarkerSize(0)

        c1.cd(i + 1)
        ROOT.gPad.SetLogx(False)

        frame = ROOT.TH2F(
            f"frame{i}",
            f"Reweighted/nominal flux ratio: {FLAVOR_LABELS[i]};Neutrino Energy (GeV);weighted / unweighted",
            100, 0, args.emax, 100, 0.0, 2.0,
        )
        frame.GetYaxis().SetTitleOffset(1.4)
        frame.Draw()

        l1 = ROOT.TLine(0, 1, args.emax, 1)
        l1.SetLineStyle(2)
        l1.SetLineColor(ROOT.kGray + 2)
        l1.Draw()

        gband.Draw("3 same")
        gcv.Draw("L same")

        leg = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.045)
        leg.AddEntry(gcv, "central value", "l")
        leg.AddEntry(gband, f"#pm1#sigma, {args.nuniv} universes", "f")
        leg.Draw()

        keep += [hnom, hcv, huniv, gcv, gband, frame, l1, leg]

    c1.cd(0)
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextSize(0.02)
    lat.DrawLatex(0.02, 0.01, f"{args.label} -- bins with zero nominal content omitted")

    c1.SaveAs(os.path.join(args.outdir, args.output))

    print("Done.")


if __name__ == "__main__":
    main()
