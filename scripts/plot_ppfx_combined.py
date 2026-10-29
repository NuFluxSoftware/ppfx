#!/usr/bin/env python3
"""Python port of plot_ppfx_combined.C (PyROOT).

Combine PPFX output across several me000z200i-style files (same universe
index = same fixed-seed systematic draw in every file, confirmed from
CentralValuesAndUncertainties::calculateParsForUniverse -- baseSeed=0, never
overridden by doReweight_dk2nu[_originalpip], so summing per-universe
histograms across files is the statistically correct combination) and remake
all three verification plots (nom-vs-cv, numu uncertainty band, ratio band)
at the combined statistics.

See notes/ppfx-weight-calculations.md and notes/ppfx-scripts-python-port.md.
"""
import argparse
import glob
import os

import ROOT

ROOT.gROOT.SetBatch(True)  # without this, batch PNG rendering subtly differs from `root -l -b -q`

from ppfx_common import COLORS, FLAVOR_LABELS, FLAVORS, style_defaults, universe_band


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--inputs", default="ppfx_out_*.root",
        help="glob pattern (or comma-separated list) of PPFX output files to combine",
    )
    ap.add_argument("--outdir", default=".", help="directory for output PNGs")
    ap.add_argument("--nuniv", type=int, default=100, help="number of PPFX universes per file")
    ap.add_argument("--pot", default=None, help="POT label for y-axis titles (default: inferred as 500,000 x nfiles)")
    ap.add_argument(
        "--label",
        default="g4numiv6 2018 production, me000z200i (ME beam, +200kA horn)",
        help="base annotation text (file count is appended automatically)",
    )
    ap.add_argument("--emax", type=float, default=60.0, help="upper edge of the ratio-plot energy frame (GeV)")
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    if "," in args.inputs:
        files = [s.strip() for s in args.inputs.split(",")]
    else:
        files = sorted(glob.glob(args.inputs))
    if not files:
        raise SystemExit(f"No input files matched: {args.inputs}")
    nfiles = len(files)
    pot = args.pot or f"{500_000 * nfiles:,}"
    label_n = f"{args.label}, {nfiles} files combined"

    style_defaults()
    fin = [ROOT.TFile(fn) for fn in files]

    # ---- sum hnom/hcv across files, per flavor ----
    hnom, hcv = [], []
    for flav in FLAVORS:
        hn = fin[0].Get(f"nom/hnom_{flav}").Clone(f"sum_hnom_{flav}")
        hc = fin[0].Get(f"nom/hcv_{flav}").Clone(f"sum_hcv_{flav}")
        hn.SetDirectory(0)
        hc.SetDirectory(0)
        for fk in fin[1:]:
            hn.Add(fk.Get(f"nom/hnom_{flav}"))
            hc.Add(fk.Get(f"nom/hcv_{flav}"))
        hnom.append(hn)
        hcv.append(hc)

    # ---- Plot 1: nominal vs CV, all 4 flavors, combined stats ----
    c1 = ROOT.TCanvas("c1", "c1", 900, 650)
    c1.SetLogy()
    leg1 = ROOT.TLegend(0.58, 0.58, 0.88, 0.88)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)

    for i in range(4):
        hnom[i].SetLineColor(COLORS[i]); hnom[i].SetLineStyle(2); hnom[i].SetLineWidth(2)
        hcv[i].SetLineColor(COLORS[i]); hcv[i].SetLineStyle(1); hcv[i].SetLineWidth(2)

    hnom[0].SetTitle(
        f"PPFX reweighting: nominal vs. central-value flux, {label_n};"
        f"Neutrino Energy (GeV);#nu / cm^{{2}} / {pot} POT / GeV"
    )
    hnom[0].GetYaxis().SetTitleOffset(1.3)
    hnom[0].SetMinimum(1e-6)
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
    lat1.SetNDC(); lat1.SetTextSize(0.03)
    lat1.DrawLatex(0.13, 0.93, f"{label_n} -- {args.nuniv} PPFX universes")
    c1.SaveAs(os.path.join(args.outdir, "ppfx_nom_vs_cv.png"))

    # ---- Plot 2: numu CV + combined many-universes band (sum over files per universe) ----
    huniv_sum = []
    for j in range(args.nuniv):
        hs = fin[0].Get(f"numu_total/htotal_numu_{j}").Clone(f"sum_htotal_numu_{j}")
        hs.SetDirectory(0)
        for fk in fin[1:]:
            hs.Add(fk.Get(f"numu_total/htotal_numu_{j}"))
        huniv_sum.append(hs)
    nbins = hcv[0].GetNbinsX()

    mean, rms = universe_band(huniv_sum, nbins)

    gband = ROOT.TGraphAsymmErrors(nbins)
    for b in range(1, nbins + 1):
        x = hcv[0].GetBinCenter(b)
        xw = hcv[0].GetBinWidth(b) / 2.0
        cvv = hcv[0].GetBinContent(b)
        lo, hi = mean[b] - rms[b], mean[b] + rms[b]
        gband.SetPoint(b - 1, x, cvv)
        gband.SetPointError(b - 1, xw, xw, max(0.0, cvv - lo), max(0.0, hi - cvv))
    gband.SetFillColorAlpha(ROOT.kBlue, 0.30)
    gband.SetLineColor(ROOT.kBlue + 1)

    c2 = ROOT.TCanvas("c2", "c2", 900, 650)
    c2.SetLogy()
    hcv[0].SetTitle(
        f"PPFX hadron-production uncertainty band: #nu_{{#mu}}, {label_n};"
        f"Neutrino Energy (GeV);#nu / cm^{{2}} / {pot} POT / GeV"
    )
    hcv[0].GetYaxis().SetTitleOffset(1.3)
    hcv[0].SetMinimum(1e-6)
    hcv[0].SetMaximum(5 * hcv[0].GetMaximum())
    hcv[0].SetLineColor(ROOT.kBlue + 2); hcv[0].SetLineWidth(2)
    hcv[0].Draw("hist")
    gband.Draw("3 same")
    hcv[0].Draw("hist same")

    leg2 = ROOT.TLegend(0.55, 0.70, 0.88, 0.88)
    leg2.SetBorderSize(0); leg2.SetFillStyle(0)
    leg2.AddEntry(hcv[0], "central value", "l")
    leg2.AddEntry(gband, f"#pm1#sigma, {args.nuniv} universes", "f")
    leg2.Draw()

    lat2 = ROOT.TLatex()
    lat2.SetNDC(); lat2.SetTextSize(0.03)
    lat2.DrawLatex(0.13, 0.93, label_n)
    c2.SaveAs(os.path.join(args.outdir, "ppfx_numu_uncertainty_band.png"))

    # ---- Plot 3: ratio vs energy, all 4 flavors, combined stats, band on ratio ----
    c3 = ROOT.TCanvas("c3", "c3", 1100, 850)
    c3.Divide(2, 2)
    keep = []

    for i, flav in enumerate(FLAVORS):
        huniv = []
        for j in range(args.nuniv):
            h = fin[0].Get(f"{flav}_total/htotal_{flav}_{j}").Clone(f"sum_ht_{flav}_{j}")
            h.SetDirectory(0)
            for fk in fin[1:]:
                h.Add(fk.Get(f"{flav}_total/htotal_{flav}_{j}"))
            huniv.append(h)

        gcv = ROOT.TGraph()
        npcv = 0
        grb = ROOT.TGraphAsymmErrors()
        npband = 0

        nb = hnom[i].GetNbinsX()
        for b in range(1, nb + 1):
            nomv = hnom[i].GetBinContent(b)
            if nomv <= 0:
                continue

            x = hnom[i].GetBinCenter(b)
            xw = hnom[i].GetBinWidth(b) / 2.0
            cvv = hcv[i].GetBinContent(b)
            gcv.SetPoint(npcv, x, cvv / nomv)
            npcv += 1

            ratios = [h.GetBinContent(b) / nomv for h in huniv]
            m = sum(ratios) / len(ratios)
            m2 = sum(r * r for r in ratios) / len(ratios)
            rms_r = max(0.0, m2 - m * m) ** 0.5

            grb.SetPoint(npband, x, m)
            grb.SetPointError(npband, xw, xw, rms_r, rms_r)
            npband += 1

        grb.SetFillColorAlpha(COLORS[i], 0.35)
        grb.SetLineColor(0)
        gcv.SetLineColor(COLORS[i]); gcv.SetLineWidth(2); gcv.SetMarkerSize(0)

        c3.cd(i + 1)
        frame = ROOT.TH2F(
            f"frame3_{i}",
            f"Reweighted/nominal flux ratio: {FLAVOR_LABELS[i]};Neutrino Energy (GeV);weighted / unweighted",
            100, 0, args.emax, 100, 0.0, 2.0,
        )
        frame.GetYaxis().SetTitleOffset(1.4)
        frame.Draw()
        l1 = ROOT.TLine(0, 1, args.emax, 1)
        l1.SetLineStyle(2); l1.SetLineColor(ROOT.kGray + 2)
        l1.Draw()
        grb.Draw("3 same")
        gcv.Draw("L same")

        leg = ROOT.TLegend(0.55, 0.75, 0.88, 0.88)
        leg.SetBorderSize(0); leg.SetFillStyle(0); leg.SetTextSize(0.045)
        leg.AddEntry(gcv, "central value", "l")
        leg.AddEntry(grb, f"#pm1#sigma, {args.nuniv} universes", "f")
        leg.Draw()

        keep += [huniv, gcv, grb, frame, l1, leg]

    c3.cd(0)
    lat3 = ROOT.TLatex()
    lat3.SetNDC(); lat3.SetTextSize(0.02)
    lat3.DrawLatex(0.02, 0.01, f"{label_n} -- bins with zero nominal content omitted")
    c3.SaveAs(os.path.join(args.outdir, "ppfx_ratio_band.png"))

    print("Done.")


if __name__ == "__main__":
    main()
