#This script is to process flux files locally and get the neutrino parent spectra. You can change the batch size, finally merge the files and while plotting take care of the normalization, as required.
#For final plot, we need to merge outputs from this and plot.
import ROOT, os, glob, re
# ROOT in batch mode (no GUI popups)
ROOT.gROOT.SetBatch(True)

# Configuration
pot_per_file = 500000.0
batch_size = 10
#input_pattern = "/pnfs/nova/persistent/stash/flux/g4numi/g4.10.04/RHC/g4numi_minervame_me000z-200i_*_0001.root"
input_pattern = "/pnfs/nova/persistent/stash/flux/g4numi/v6r3/old_target/me000z200i/g4numiv6_minervame_me000z200i_*_0200.root"


# Sort files by embedded number
def extract_number(filename):
    match = re.search(r"_(\d+)_0200\.root", filename)
    return int(match.group(1)) if match else -1

fluxpaths = sorted(glob.glob(input_pattern), key=extract_number)
print("Total files found: {}".format(len(fluxpaths)))

# Process files in batches
for batch_start in range(0, len(fluxpaths), batch_size):
    batch = fluxpaths[batch_start:batch_start + batch_size]
    print("\nProcessing batch {} to {}".format(batch_start + 1, batch_start + len(batch)))

    fluxchain = ROOT.TChain("dk2nuTree")
    for file in batch:
        #print("Adding file {}:".format(file))
        fluxchain.Add(file)

    # POT Calculation
    #total_pot = pot_per_file * len(batch)

    # Histograms
    h_muon_par     = ROOT.TH1F("h_muon_par{}".format(batch_start),     "#nu from #mu^{-}",        320, 0, 20)
    h_kaon_par = ROOT.TH1F("h_kaon_par{}".format(batch_start), "#nu from K^{+}",   320, 0, 20)
  h_kl_par      = ROOT.TH1F("h_kl_par{}".format(batch_start),      "#nu from K^{0}",           320, 0, 20)
  h_pion_par  = ROOT.TH1F("h_pion_par{}".format(batch_start),  "#nu from pi^{+}",     320, 0, 20)

    histogram = [h_muon_par, h_kaon_par, h_kl_par, h_pion_par]
    colors = [ROOT.kBlack, ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]

    for i, hist in enumerate(histogram):
        hist.SetLineColor(colors[i])
        hist.SetLineWidth(2)
        hist.GetXaxis().SetTitle("Neutrino Energy (GeV)")
        hist.GetYaxis().SetTitle("Neutrinos/m^{2}/10^{8}POT")
        hist.GetXaxis().CenterTitle(True)
        hist.GetYaxis().CenterTitle(True)

    # Fill histograms (adjust variable names as per your ROOT files)
    fluxchain.Draw("dk2nu.nuray.E[3]>>h_muon_par{}".format(batch_start),"dk2nu.decay.nimpwt*dk2nu.nuray.wgt[3]/TMath::Pi()*(dk2nu.decay.ptype==-13)")
    fluxchain.Draw("dk2nu.nuray.E[3]>>h_kaon_par{}".format(batch_start),"dk2nu.decay.nimpwt*dk2nu.nuray.wgt[3]/TMath::Pi()*(dk2nu.decay.ptype==321)")
    fluxchain.Draw("dk2nu.nuray.E[3]>>h_kl_par{}".format(batch_start),"dk2nu.decay.nimpwt*dk2nu.nuray.wgt[3]/TMath::Pi()*(dk2nu.decay.ptype==130 || dk2nu.decay.ptype==310)")
    fluxchain.Draw("dk2nu.nuray.E[3]>>h_pion_par{}".format(batch_start),"dk2nu.decay.nimpwt*dk2nu.nuray.wgt[3]/TMath::Pi()*(dk2nu.decay.ptype==211)")
    # Scale histograms
   # for hist in histogram:
   #     hist.Scale(1.0 / total_pot, "width")

    # Output root file
    output_root = "output_batch_{}_to_{}.root".format(batch_start + 1, batch_start + len(batch))
    f_out = ROOT.TFile(output_root, "RECREATE")
    ROOT.gStyle.SetOptStat(0)

    for hist in histogram:
        hist.Write()

  
    canvas = ROOT.TCanvas("canvas", "Neutrino parent  Spectrum", 800, 600)
    h_muon_par.Draw("HIST")
    h_kaon_mu_par.Draw("HIST SAME")
    h_kl_par.Draw("HIST SAME")
    h_pion_par.Draw("HIST SAME")

    legend = ROOT.TLegend(0.6, 0.7, 0.88, 0.88)
    legend.AddEntry(h_muon_par, "{#mu}", "l")
    legend.AddEntry(h_kaon_par, "Kaon", "l")
    legend.AddEntry(h_kl_par, "K_{L}", "l")
    legend.AddEntry(h_pion_par, "#pi", "l")
    legend.Draw()

    text = ROOT.TLatex()
    text.SetNDC()
    text.SetTextSize(0.04)
    text.SetTextColor(ROOT.kGray + 2)
   # text.DrawLatex(0.7, 0.93, "NOvA Simulation")

    canvas.Write("combined_canvas")
    canvas.SaveAs(output_root.replace(".root", ".png"))

    
    f_out.Close()

print("\nAll batches processed successfully.")
