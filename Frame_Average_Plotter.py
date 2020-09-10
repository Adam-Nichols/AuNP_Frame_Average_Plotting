import numpy as np
import pandas as pd
import math
import glob
import matplotlib.pyplot as plt
from os import chdir

# Settings #
plot_radius = int(25)
plot_atoms = ['Au', 'Water', 'PEG']  # selection of atoms below
# Full atom list: ['Cr', 'NR', 'OR', 'SP', 'HCM', 'HNR', 'OT', 'HT', 'Au', 'Water', 'PEG']
# Individual plot type list for fast copy/paste: ['Cr', 'NR', 'OR', 'SP', 'Au', 'Water']
# Combined plot type list for fast copy/paste: ['Au', 'Water', 'PEG']
chdir(r"C:\MBN\Au135_08PEG2_Water_equil_frames_test")


read_files = glob.glob("*.xyz")
file_count = len(read_files)

complete_data = []

for i in range(file_count):
    data_sheet = open(read_files[i], "r")
    data_sheet = data_sheet.read().split()
    data_sheet = data_sheet[4:]
    complete_data.append(data_sheet)

raw_positions = []
for sublist in complete_data:
    for item in sublist:
        raw_positions.append(item)

au_count = int(raw_positions.count("Au")*(1/file_count))

y_labels = []
radius = []
data_for_frames = []

i = 0
while i < len(raw_positions):
    data_for_frames.append(raw_positions[i:i+4])
    i += 4

for i in range(int(len(raw_positions)/4)):
    y_labels.append(raw_positions[i*4])

data = pd.DataFrame(data_for_frames, columns=["Type", "x", "y", "z"])
data.set_index("Type", inplace=True)
data = data.astype({"x": float, "y": float, "z": float})

atom_types = list(data.index.unique())

row_count = data.shape[0]
row_count = int(row_count*(1/file_count))


for i in range(file_count):
    au_data = data.iloc[(i+1)*row_count-au_count:(i+1)*row_count, :]
    x_avg = au_data.mean()[0]
    y_avg = au_data.mean()[1]
    z_avg = au_data.mean()[2]
    data["x"][i*row_count:i*row_count+row_count] -= x_avg
    data["y"][i*row_count:i*row_count+row_count] -= y_avg
    data["z"][i*row_count:i*row_count+row_count] -= z_avg

data["Distance"] = np.sqrt(data["x"]**2 + data["y"]**2 + data["z"]**2)
data["Shell"] = data["Distance"].apply(np.floor)
data = data.astype({"Shell": int})

density = pd.DataFrame()
plot_radius += 1
density["Radius"] = range(plot_radius)
density["Volume"] = ((4/3)*math.pi*(density["Radius"]**3))-((4/3)*math.pi*((density["Radius"]-1)**3))
density.set_index("Radius", inplace=True)

for f in atom_types:
    density[f] = 0

for i in range(plot_radius):
    for f in atom_types:
        try:
            series = data.loc[[f], ["Shell"]] == i
            count = series.Shell.value_counts()[1]
            density.at[i+1, f] = count
        except KeyError:
            pass

for f in atom_types:
    density[f] = density[f]*(1/file_count)

density["Water"] = (density["OT"] + density["HT"])*(1/3)
atom_types.append("Water")
density["PEG"] = density["Cr"] + density["NR"] + density["OR"] + density["SP"] + density["HCM"] + density["HNR"]
atom_types.append("PEG")
density.loc[2, "Au"] += density.loc[1, "Au"]
density.loc[1, "Au"] = 0
density = density[:-1]
density["Au"] = density["Au"]*(1/3)

for f in atom_types:
    density[f+"_rho"] = density[f]/density["Volume"]

x_data = range(len(density))
for f in plot_atoms:
    y_data = density[f+"_rho"]
    if f == "Au":
        plt.plot(x_data, y_data, "-o", label=f+" (1/3)")
    else:
        plt.plot(x_data, y_data, "-o", label=f)
    plt.legend()

plt.xlabel("Radius (Angstroms)")
plt.ylabel("Number Density")
plt.hlines(0.0334, 0, len(density), colors="gray", linestyles="dashed")
plt.show()

density.rename(columns={"Au": "Au_(1/3)", "Au_rho": "Au_(1/3)_rho"})
density.to_csv("Frame_Averaged_Densities.csv")