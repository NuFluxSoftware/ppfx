"""Shared helpers for the PPFX interaction-map (IMap) plotting scripts
(plot_flux_map.py, hp_coverage.py).

get_h1/get_h2 are the only two things these scripts actually used from the
external ~/root_utils/move_stats.C utility library -- reimplemented here
(open a file, fetch a histogram, detach it from the file) rather than
vendoring/porting that whole library, per user request. `center_titles()`
isn't factored out the same way: each of its call sites only ever acts on
the one histogram just drawn in the current pad, so the callers just call
`h.GetXaxis().CenterTitle()`/`h.GetYaxis().CenterTitle()` directly instead
of going through a generic gPad-primitives-iterating utility.
"""
import ROOT


def get_h1(hname, fname):
    f = ROOT.TFile(fname)
    h = f.Get(hname)
    h.SetDirectory(0)
    return h


def get_h2(hname, fname):
    return get_h1(hname, fname)  # TFile::Get + SetDirectory(0) is identical for TH1/TH2
