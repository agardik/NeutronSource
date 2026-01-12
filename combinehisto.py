import ROOT

# Prevent ROOT from taking ownership
ROOT.TH1.AddDirectory(False)

# Open ROOT file
f = ROOT.TFile.Open("/home/kali/Software/Newtron source git/NeutronSource/Results/build/DoseResults_EnergyDep.root", "READ")
if not f or f.IsZombie():
    raise RuntimeError("Cannot open DoseResults_EnergyDep.root")

# Get histograms
h_alpha = f.Get("alphadep")
h_li    = f.Get("Li7dep")
h_all   = f.Get("all_dep")

if not h_alpha or not h_li or not h_all:
    raise RuntimeError("One or more histograms not found")

# Normalize for fair comparison (optional but recommended)
#for h in (h_alpha, h_li, h_all):
    #if h.Integral() > 0:
        #h.Scale(1.0 / h.Integral())

max_val = max(h_alpha.GetMaximum(), h_li.GetMaximum(), h_all.GetMaximum())
h_alpha.SetMaximum(max_val * 0.1)  # 20% extra space above the highest peak
h_alpha.SetMinimum(0)


# Style
h_alpha.SetTitle("Energy Deposition")
h_alpha.SetLineColor(ROOT.kRed)
h_li.SetLineColor(ROOT.kBlue)
h_all.SetLineColor(ROOT.kBlack)

for h in (h_alpha, h_li, h_all):
    h.SetLineWidth(2)

# Canvas
c = ROOT.TCanvas("c", "Energy Deposition Overlay", 900, 700)

# Axis labels
h_alpha.GetXaxis().SetTitle("Energy deposited [MeV]")
h_alpha.GetYaxis().SetTitle("Counts")

# Disable stats box
h_alpha.SetStats(False)
h_li.SetStats(False)
h_all.SetStats(False)

# Draw
h_alpha.Draw("HIST")
h_li.Draw("HIST SAME")
h_all.Draw("HIST SAME")

# Legend
leg = ROOT.TLegend(0.65, 0.7, 0.88, 0.88)
leg.AddEntry(h_alpha, "Alpha", "l")
leg.AddEntry(h_li, "Li-7", "l")
leg.AddEntry(h_all, "All particles", "l")
leg.Draw()

# Optional: log scale (useful for tails)
# c.SetLogy()

# Save
c.SaveAs("EnergyDep_overlay.png")
input("Press Enter to exit...")
f.Close()
