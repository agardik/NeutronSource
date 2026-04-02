import ROOT
import math
from collections import defaultdict
import matplotlib.pyplot as plt

# Open ROOT file
file = ROOT.TFile.Open("/media/sf_vm-share/Neutron source output/rootoutput/testoutput.root")

# Get tree
tree = file.Get("HeliumEdep")

# Strip width
W = 32 / 1  # change to desired width

# Threshold
threshold = 1e-9

# Histogram parameters
xmin = 20
xmax = 180
nbins = int((xmax - xmin) / W)

# Map: (event, strip) -> accumulated Edep
stripEventEdep = defaultdict(float)

# Loop over entries
nentries = tree.GetEntries()
for i in range(nentries):
    tree.GetEntry(i)
    fEvent = getattr(tree, "fEvent")
    Edep = getattr(tree, "Edep")
    fX = getattr(tree, "fX")

    if Edep < threshold:
        continue

    # Determine strip index
    stripID = math.floor((fX - xmin) / W)
    key = (fEvent, stripID)
    stripEventEdep[key] += Edep

# Fill histogram data
x_positions = []
y_edep = []

for (event, stripID), edep in stripEventEdep.items():
    if edep < threshold:
        continue
    xpos = xmin + (stripID + 0.5) * W  # center of strip
    x_positions.append(xpos)
    y_edep.append(edep)

# Normalize histogram
total_edep = sum(y_edep)
if total_edep > 0:
    y_edep = [y / total_edep for y in y_edep]

# Plot using matplotlib
plt.figure(figsize=(10, 6))
plt.bar(x_positions, y_edep, width=W*0.9, align='center', edgecolor='black')
plt.xlabel("X position")
plt.ylabel("Normalized Edep")
plt.title("Energy deposition per strip")

# Optional: Add top axis for Strip ID
ax = plt.gca()
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
nStrips = nbins
strip_labels = range(nStrips)
strip_positions = [xmin + (i + 0.5) * W for i in strip_labels]
ax2.set_xticks(strip_positions)
ax2.set_xticklabels(strip_labels)
ax2.set_xlabel("Strip ID")

plt.tight_layout()
plt.savefig("strip_edep.png")
plt.show()

print(f"Total strip-event combinations: {len(stripEventEdep)}")
print(f"Histogram integral = {sum(y_edep)}")
print(f"Sum of weights = {sum(y_edep)}")