from rdkit import Chem
import tmap as tm
from map4 import MAP4Calculator

# Function to read SMILES from a file and categorize them
def read_smiles(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return [line.strip() for line in lines]

# Initialize MAP4 and Minhash
dim = 1024
MAP4 = MAP4Calculator(dimensions=dim)
ENC = tm.Minhash(dim)

# Read SMILES from file
smiles_list = read_smiles("smiles.smi")

# Split the list into pep and str
peps = smiles_list[:12]  # First 12 are pep1-12
strs = smiles_list[12:]  # Remaining are str1-9

# Function to calculate MAP4 fingerprint
def calculate_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return MAP4.calculate(mol)

# Calculate fingerprints for each molecule
pep_fps = [calculate_fingerprint(smiles) for smiles in peps]
str_fps = [calculate_fingerprint(smiles) for smiles in strs]

with open("minhash_distance.txt", "w") as file:
    # Compute distances between str and pep
    for i, str_fp in enumerate(str_fps):
        for j, pep_fp in enumerate(pep_fps):
            distance = 1 - ENC.get_distance(str_fp, pep_fp) # here we transform the dissimilaty into similarity
            file.write(f"str{i+1} to pep{j+1}: {distance}\n")

    # Compute distances among str (excluding self comparisons)
    for i, str_fp1 in enumerate(str_fps):
        for j, str_fp2 in enumerate(str_fps):
            if i != j:
                distance = 1 - ENC.get_distance(str_fp1, str_fp2)
                file.write(f"str{i+1} to str{j+1}: {distance}\n")


import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

# Read and process the data from the file
def read_data(file_path):
    str_to_pep = {f"str{i}": [] for i in range(1, 10)}
    str_to_str = {f"str{i}": [] for i in range(1, 10)}

    with open(file_path, "r") as file:
        for line in file:
            parts = line.strip().split(": ")
            pair, distance = parts[0], float(parts[1])
            str_molecule, other_molecule = pair.split(" to ")
            if other_molecule.startswith("pep"):
                str_to_pep[str_molecule].append(distance)
            else:
                str_to_str[str_molecule].append(distance)

    return str_to_pep, str_to_str

# Calculate averages and standard deviations
def calculate_stats(data):
    averages = []
    std_devs = []
    for key in sorted(data.keys()):
        values = data[key]
        averages.append(np.mean(values))
        std_devs.append(np.std(values))
    return averages, std_devs

# File path
file_path = "minhash_distance.txt"

# Process the data
str_to_pep, str_to_str = read_data(file_path)

# Calculate stats
pep_averages, pep_std_devs = calculate_stats(str_to_pep)
str_averages, str_std_devs = calculate_stats(str_to_str)

# Create bar charts
x = np.arange(len(str_to_pep))  # the label locations

fig, ax = plt.subplots()
bar_width = 0.35  # the width of the bars

# Set global font properties
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 12  # You can adjust the size as needed


fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the size of the figure

rects1 = ax.bar(x - bar_width/2, pep_averages, bar_width, yerr=pep_std_devs, label='Average similarity to Pep', capsize=2, color='skyblue')
rects2 = ax.bar(x + bar_width/2, str_averages, bar_width, yerr=str_std_devs, label='Average similarity to Str', capsize=2, color='salmon')

ax.set_ylabel('Average Jaccard Similarity')
ax.set_xticks(x)
ax.set_xticklabels([f"str{i}" for i in range(1, 10)], rotation=45, ha="right")  # Rotate labels and adjust alignment

ax.legend()
plt.tight_layout()  # Adjust layout

plt.savefig("minhash_distance_similarity.svg", format='svg')
plt.show()