"""Shared constants/helpers for the PPFX plotting scripts (plot_ppfx_weights.py,
plot_ppfx_ratio.py, plot_ppfx_combined.py). Python ports of plot_ppfx_*.C,
see notes/ppfx-scripts-python-port.md.
"""
import ROOT

FLAVORS = ["numu", "numubar", "nue", "nuebar"]
FLAVOR_LABELS = ["#nu_{#mu}", "#bar{#nu}_{#mu}", "#nu_{e}", "#bar{#nu}_{e}"]
COLORS = [ROOT.kBlue + 1, ROOT.kRed + 1, ROOT.kGreen + 2, ROOT.kMagenta + 1]


def style_defaults():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadBottomMargin(0.12)


def universe_band(universe_hists, nbins, divisor_hist=None):
    """Bin-by-bin mean/RMS across a list of per-universe TH1s.

    If divisor_hist is given, each universe's bin content is first divided by
    the divisor histogram's bin content at that bin (the ratio-band case);
    otherwise the raw bin content is used (the absolute-band case). Returns
    (mean, rms) arrays, 1-indexed by bin (index 0 unused) to match ROOT's
    GetBinContent(b) convention.
    """
    nuniv = len(universe_hists)
    mean = [0.0] * (nbins + 1)
    rms = [0.0] * (nbins + 1)
    for b in range(1, nbins + 1):
        div = divisor_hist.GetBinContent(b) if divisor_hist is not None else 1.0
        s, s2 = 0.0, 0.0
        for h in universe_hists:
            v = h.GetBinContent(b) / div if divisor_hist is not None else h.GetBinContent(b)
            s += v
            s2 += v * v
        m = s / nuniv
        mean[b] = m
        rms[b] = (max(0.0, s2 / nuniv - m * m)) ** 0.5
    return mean, rms
