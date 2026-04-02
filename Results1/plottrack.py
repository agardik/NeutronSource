import ROOT
from collections import defaultdict

# -------------------------
# USER SETTINGS
# -------------------------
root_file = "/home/kali/Software/Newtron source git/NeutronSource/Results/build/testoutput.root"
tree_name = "HeliumEdep"   # CHANGE if needed

E_MIN = 1.3   # MeV
E_MAX = 1.5   # MeV

MAX_TRACKS = 20   # limit for readability
PLANE = "xz"      # "xz" or "xy" or "yz"

# -------------------------
# OPEN FILE
# -------------------------
f = ROOT.TFile.Open(root_file)
tree = f.Get(tree_name)
tree.GetEntry(0)  # to load branches
tree.Show(0)
# -------------------------
# FIRST PASS: total Edep per event (alphas only)
# -------------------------
event_edep = defaultdict(float)

for entry in tree:
    if entry.fParticleName != "alpha":
        continue
    event_edep[entry.fEvent] += entry.Edep

# Select events in energy band
selected_events = {
    evt for evt, E in event_edep.items()
    if E_MIN <= E <= E_MAX
}

print(f"Selected {len(selected_events)} alpha events in energy band")

# -------------------------
# SECOND PASS: collect track points
# -------------------------
tracks = defaultdict(list)

for entry in tree:
    if entry.fParticleName != "alpha":
        continue
    if entry.fEvent not in selected_events:
        continue

    tracks[entry.fEvent].append((entry.fX, entry.fY, entry.fZ))

# Limit number of tracks
tracks = dict(list(tracks.items())[:MAX_TRACKS])

# -------------------------
# PLOTTING
# -------------------------
canvas = ROOT.TCanvas("c", "Alpha Tracks", 900, 700)
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)

graphs = []
colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2,
          ROOT.kMagenta, ROOT.kCyan+2, ROOT.kOrange+7]

for i, (evt, points) in enumerate(tracks.items()):
    xs, ys = [], []

    for x, y, z in points:
        if PLANE == "xz":
            xs.append(z)
            ys.append(x)
        elif PLANE == "xy":
            xs.append(x)
            ys.append(y)
        elif PLANE == "yz":
            xs.append(z)
            ys.append(y)

    gr = ROOT.TGraph(len(xs))
    for j in range(len(xs)):
        gr.SetPoint(j, xs[j], ys[j])

    gr.SetLineColor(colors[i % len(colors)])
    gr.SetLineWidth(2)
    gr.SetTitle("")

    if i == 0:
        gr.Draw("AL")
        gr.GetXaxis().SetTitle(PLANE[0] + " (cm)")
        gr.GetYaxis().SetTitle(PLANE[1] + " (cm)")
    else:
        gr.Draw("L SAME")

    legend.AddEntry(gr, f"Event {evt}", "l")
    graphs.append(gr)

legend.Draw()

canvas.SaveAs("alpha_tracks_energy_band.png")
canvas.SaveAs("alpha_tracks_energy_band.root")
