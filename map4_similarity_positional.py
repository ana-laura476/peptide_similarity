import pandas as pd
from rdkit import Chem
from map4 import MAP4Calculator
import tmap as tm
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

# Set global font properties
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 12  # You can adjust the size as needed

# Load the CSV file
df = pd.read_csv("converted_sequences.csv", sep = ';')
print(df)

# Initialize MAP4 and Minhash
dim = 1024
map4_calc = MAP4Calculator(dimensions=dim)
mh = tm.Minhash(dim)

# Function to calculate MAP4 fingerprint
def calculate_map4_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return map4_calc.calculate(mol)

# Process the data
results = []

# Iterate through each 'str' entry
for str_index in range(1, 10):  # str1 to str9
    str_name = f"str{str_index}"
    for position in range(12):  # 12 positions in the peptide
        # Get the SMILES string for the current position from the current 'str'
        ref_smiles = eval(df.loc[df['name'] == str_name, 'smiles'].values[0])[position]
        ref_fp = calculate_map4_fingerprint(ref_smiles)
        
        # Compare with all other peptides' corresponding position
        for index, row in df.iterrows():
            if row['name'] != str_name:  # Exclude self-comparison
                target_smiles = eval(row['smiles'])[position]
                target_fp = calculate_map4_fingerprint(target_smiles)
                distance = mh.get_distance(ref_fp, target_fp)
                similarity = distance  # Convert distance to similarity
                
                # Store the results
                results.append({
                    'reference': str_name,
                    'target': row['name'],
                    'position': position + 1,
                    'minhash distance': similarity
                })

# Convert results to a DataFrame
results_df = pd.DataFrame(results)

# Save results to a new CSV file
results_df.to_csv("amino_acid_position_similarity.csv", index=False)

print("Similarity calculations complete. Results saved to amino_acid_position_similarity.csv")

# Load the results CSV file
results_df = pd.read_csv("amino_acid_position_similarity.csv")
fig = plt.figure(figsize=(20, 40))  # Adjust the size for better fit if necessary
outer_grid = gridspec.GridSpec(9, 1, figure=fig, hspace=0.8)

# Iterate over each 'str' and plot in a composite subplot
for str_index in range(1, 10):
    str_name = f"str{str_index}"
    str_results = results_df[results_df['reference'] == str_name]
    positions = range(1, 13)  # 12 positions

    # Create a subplot for 'pep' above and 'str' below
    inner_grid = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer_grid[str_index-1], hspace=0.3)
    ax_pep = fig.add_subplot(inner_grid[0])
    ax_str = fig.add_subplot(inner_grid[1], sharex=ax_pep)

    # Calculate averages and standard deviations
    avg_similarity_pep = []
    std_dev_pep = []
    avg_similarity_str = []
    std_dev_str = []

    for pos in positions:
        pep_similarities = 1- (str_results[str_results['target'].str.startswith('pep') & (str_results['position'] == pos)]['minhash distance'])
        avg_similarity_pep.append(pep_similarities.mean())
        std_dev_pep.append(pep_similarities.std())

        str_similarities = 1- (str_results[str_results['target'].str.startswith('str') & (str_results['target'] != str_name) & (str_results['position'] == pos)]['minhash distance'])
        avg_similarity_str.append(str_similarities.mean())
        std_dev_str.append(str_similarities.std())

    # Plot for 'pep'
    x = np.arange(len(positions))  # the label locations
    width = 0.35  # the width of the bars
    bars_pep = ax_pep.bar(x, avg_similarity_pep, width, yerr=std_dev_pep, color='skyblue', capsize=2)
    bars_str = ax_str.bar(x, avg_similarity_str, width, yerr=std_dev_str, color='salmon', capsize=2)
    ax_pep.bar(x, avg_similarity_pep, width, yerr=std_dev_pep, color='skyblue')
    ax_pep.set_ylim([0, 1])  # Fixed y-axis limits

    # Plot for 'str'
    ax_str.bar(x, avg_similarity_str, width, yerr=std_dev_str, color='salmon')
    ax_str.set_ylim([0, 1])  # Fixed y-axis limits

    # Set labels and titles
    ax_pep.set_title(f'{str_name} Similarity', fontsize=14)
    ax_str.set_xlabel('Amino acid position', fontsize=12)
    # Add a shared legend to the right of the subplots
    handles_pep, labels_pep = ax_pep.get_legend_handles_labels()
    handles_str, labels_str = ax_str.get_legend_handles_labels()
    fig.legend(handles_pep + handles_str, ['Average similarity to Pep', 'Average similarity to Str'], 
               loc='upper right', bbox_to_anchor=(1.1, 1.25), fontsize=12, frameon=False)

    # Add a shared y-axis label
    ax_label = fig.add_subplot(outer_grid[str_index-1], frame_on=False)  # Add a new axis for the label
    ax_label.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)  # Hide tick marks
    ax_label.set_ylabel('Average Jaccard Similarity', labelpad=20)  # Set the y-axis label
    ax_label.yaxis.set_label_position("left")  # Position the label to the left
    ax_label.yaxis.set_ticks_position('none') 

    ax_str.set_xticks(x)
    ax_str.set_xticklabels([f'{i}' for i in positions], fontsize=10)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    plt.setp(ax_pep.get_xticklabels(), visible=False)

    # Add a common y-axis label
    fig.text(0.04, 0.5, 'Similarity', va='center', rotation='vertical', fontsize=14)

# Add a single legend for the entire figure

# Adjust the layout
fig.tight_layout(rect=[0.05, 0.03, 1, 0.95])

# Save the figure as SVG
plt.savefig("minhash_distance_similarity_12aminoacids.svg", format='svg')
plt.show()