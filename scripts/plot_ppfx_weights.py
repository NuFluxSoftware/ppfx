#!/usr/bin/env python3
"""Python port of plot_ppfx_weights.C (PyROOT).

Plot 1: nominal vs. central-value-reweighted flux, all 4 flavors.
Plot 2: numu central value + many-universes +/-1sigma uncertainty band.

See notes/ppfx-weight-calculations.md and notes/ppfx-scripts-python-port.md.
"""
import argparse
import os

import ROOT

ROOT.gROOT.SetBatch(True)  # without this, batch PNG rendering subtly differs from `root -l -b -q`

from ppfx_common import COLORS, FLAVOR_LABELS, FLAVORS, style_defaults, universe_band


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", default="ppfx_out_0000.root", help="PPFX doReweight_dk2nu output file")
    ap.add_argument("--outdir", default=".", help="directory for output PNGs")
    ap.add_argument("--nuniv", type=int, default=100, help="number of PPFX universes in the file")
    ap.add_argument("--pot", default="500,000", help="POT label for the y-axis title")
    ap.add_argument(
        "--label",
        default="g4numiv6 2018 production, me000z200i (ME beam, +200kA horn), file 0000",
        help="top-left annotation text (universe count is appended automatically)",
    )
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    style_defaults()
    f = ROOT.TFile(args.input)

    # ---- Plot 1: nominal vs. central-value-reweighted flux, all 4 flavors ----
    c1 = ROOT.TCanvas("c1", "c1", 900, 650)
    c1.SetLogy()

    leg1 = ROOT.TLegend(0.58, 0.58, 0.88, 0.88)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)

    hnom, hcv = [], []
    for i, flav in enumerate(FLAVORS):
        hn = f.Get(f"nom/hnom_{flav}")
        hc = f.Get(f"nom/hcv_{flav}")
        hn.SetDirectory(0)
        hc.SetDirectory(0)
        hn.SetLineColor(COLORS[i])
        hn.SetLineStyle(2)
        hn.SetLineWidth(2)
        hc.SetLineColor(COLORS[i])
        hc.SetLineStyle(1)
        hc.SetLineWidth(2)
        hnom.append(hn)
        hcv.append(hc)

    hnom[0].SetTitle(
        f"PPFX reweighting: nominal vs. central-value flux, {args.label};"
        f"Neutrino Energy (GeV);#nu / cm^{{2}} / {args.pot} POT / GeV"
    )
    hnom[0].GetYaxis().SetTitleOffset(1.3)
    hnom[0].SetMinimum(1e-7)
    hnom[0].SetMaximum(5 * hnom[0].GetMaximum())
    hnom[0].Draw("hist")
    hcv[0].Draw("hist same")
    for i in range(1, 4):
        hnom[i].Draw("hist same")
        hcv[i].Draw("hist same")

    for i in range(4):
        leg1.AddEntry(hcv[i], f"{FLAVOR_LABELS[i]} (CV reweighted)", "l")
    leg1.AddEntry(hnom[0], "nominal (dashed, any flavor)", "l")
    leg1.Draw()

    lat1 = ROOT.TLatex()
    lat1.SetNDC()
    lat1.SetTextSize(0.03)
    lat1.DrawLatex(0.13, 0.93, f"{args.label} -- {args.nuniv} PPFX universes available")
    c1.SaveAs(os.path.join(args.outdir, "ppfx_nom_vs_cv.png"))

    # ---- Plot 2: numu central value + many-universes 1-sigma uncertainty band ----
    huniv = [f.Get(f"numu_total/htotal_numu_{j}") for j in range(args.nuniv)]
    nbins = huniv[0].GetNbinsX()

    mean, rms = universe_band(huniv, nbins)

    gband = ROOT.TGraphAsymmErrors(nbins)
    for b in range(1, nbins + 1):
        x = hcv[0].GetBinCenter(b)
        xw = hcv[0].GetBinWidth(b) / 2.0
        lo, hi = mean[b] - rms[b], mean[b] + rms[b]
        cvv = hcv[0].GetBinContent(b)
        gband.SetPoint(b - 1, x, cvv)
        gband.SetPointError(b - 1, xw, xw, max(0.0, cvv - lo), max(0.0, hi - cvv))
    gband.SetFillColorAlpha(ROOT.kBlue, 0.30)
    gband.SetLineColor(ROOT.kBlue + 1)

    c2 = ROOT.TCanvas("c2", "c2", 900, 650)
    c2.SetLogy()

    hcv[0].SetTitle(
        f"PPFX hadron-production uncertainty band: #nu_{{#mu}}, {args.label};"
        f"Neutrino Energy (GeV);#nu / cm^{{2}} / {args.pot} POT / GeV"
    )
    hcv[0].GetYaxis().SetTitleOffset(1.3)
    hcv[0].SetMinimum(1e-7)
    hcv[0].SetMaximum(5 * hcv[0].GetMaximum())
    hcv[0].SetLineColor(ROOT.kBlue + 2)
    hcv[0].SetLineWidth(2)
    hcv[0].Draw("hist")
    gband.Draw("3 same")
    hcv[0].Draw("hist same")  # redraw CV line on top of band

    leg2 = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.AddEntry(hcv[0], "central value", "l")
    leg2.AddEntry(gband, f"#pm1#sigma, {args.nuniv} reweighted universes", "f")
    leg2.Draw()

    lat2 = ROOT.TLatex()
    lat2.SetNDC()
    lat2.SetTextSize(0.03)
    lat2.DrawLatex(0.13, 0.93, args.label)
    c2.SaveAs(os.path.join(args.outdir, "ppfx_numu_uncertainty_band.png"))

    print("Done.")


if __name__ == "__main__":
    main()
