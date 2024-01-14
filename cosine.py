import sklearn
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


# Convert the texts into TF-IDF vectors
pep1 = 'cysteine tryptophan leucine asparagine histidine proline glutamine glycine proline proline serine cysteine'
pep2 = 'cysteine tryptophan leucine asparagine histidine proline glutamine alanine proline proline glutamate cysteine'
pep3 = 'cysteine histidine proline glutamine glycine aspartate arginine tyrosine glutamate lysine glutamate cysteine'
pep4 = 'cysteine serine leucine threonine asparagine histidine proline glutamine asparagine phenylalanine asparagine cysteine'
pep5 = 'cysteine isoleucine glutamate leucine tryptophan arginine histidine proline glutamine glycine proline cysteine'
pep6 = 'cysteine histidine proline glutamine aspartate asparagine arginine arginine threonine phenylalanine glycine cysteine'
pep7 = 'cysteine histidine proline glutamine aspartate asparagine arginine asparagine asparagine serine serine cysteine'
pep8 = 'cysteine histidine proline glutamine aspartate asparagine arginine alanine asparagine serine serine cysteine'
pep9 = 'cysteine isoleucine histidine methionine phenylalanine asparagine histidine proline glutamine asparagine aspartate cysteine'
pep10 = 'cysteine isoleucine proline asparagine glutamate aspartate histidine proline glutamine asparagine serine cysteine'
pep11 = 'cysteine isoleucine histidine methionine phenylalanine asparagine histidine proline glutamine asparagine aspartate cysteine'
pep12 = 'cysteine proline histidine leucine serine lysine alanine histidine proline glutamine valine cysteine'

str1 = 'cysteine proline isoleucine asparagine isoleucine phenylalanine histidine proline proline proline aspartate cysteine'
str2 = 'cysteine proline phenylalanine asparagine isoleucine lysine asparagine proline aspartate aspartate glycine cysteine'
str3 = 'cysteine proline tyrosine asparagine isoleucine alanine glycine proline proline leucine tyrosine cysteine'
str4 = 'cysteine proline phenylalanine asparagine isoleucine leucine valine proline phenylalanine isoleucine aspartate cysteine'
str5 = 'cysteine proline tyrosine asparagine isoleucine isoleucine asparagine proline threonine phenylalanine aspartate cysteine'
str6 = 'cysteine proline phenylalanine asparagine isoleucine alanine serine proline phenylalanine tyrosine glycine cysteine'
str7 = 'cysteine proline isoleucine asparagine valine glutamine tyrosine serine threonine glycine threonine cysteine'
str8 = 'cysteine proline tyrosine asparagine valine glutamate arginine proline aspartate serine histidine cysteine'
str9 = 'cysteine proline phenylalanine asparagine phenylalanine phenylalanine isoleucine tryptophan asparagine aspartate glutamate cysteine'


# Convert the texts into TF-IDF vectors
vectorizer = TfidfVectorizer()
vectors = vectorizer.fit_transform([pep1, pep2, pep3, pep4, pep5, pep6, pep7, pep8, pep9, pep10, pep11, pep12,
str1, str2, str3, str4, str5, str6, str7, str8, str9])

cosine_similarities  = cosine_similarity(vectors)
print(cosine_similarities )

# Open a new text file for writing
with open('cosine_similarities.txt', 'w') as file:
    # Iterate over each row in the cosine similarity matrix
    for row in cosine_similarities:
        # Format each number to have only three decimal points and write the row to the file
        file.write(','.join(['{:.3f}'.format(value) for value in row]) + '\n')


# Calculate the average cosine similarities for Str sequences with both Pep and Str sequences
# Initialize lists to store the average similarities and standard deviations
avg_similarities_pep = []
std_dev_pep = []
avg_similarities_str = []
std_dev_str = []

# Loop through each Str sequence
for i in range(9):  # There are 9 Str sequences
    str_index = i + 12  # Indices 12-20 in the cosine similarity matrix correspond to Str1-9

    # Similarities with Pep sequences (first 12)
    similarities_with_pep = cosine_similarities[str_index, :12]
    avg_similarities_pep.append(np.mean(similarities_with_pep))
    std_dev_pep.append(np.std(similarities_with_pep))

    # Similarities with Str sequences (excluding self-similarity)
    similarities_with_str = np.concatenate([cosine_similarities[str_index, 12:str_index], 
                                            cosine_similarities[str_index, str_index+1:]])
    avg_similarities_str.append(np.mean(similarities_with_str))
    std_dev_str.append(np.std(similarities_with_str))

# Creating the plot
str_labels = [f'str{i+1}' for i in range(9)]
x = np.arange(len(str_labels))  # the label locations
width = 0.35  # the width of the bars
fig, ax = plt.subplots(figsize=(10, 6))  # Adjust the size of the figure

rects1 = ax.bar(x - width/2, avg_similarities_pep, width, yerr=std_dev_pep, label='Average similarity to Pep', capsize=2, color='skyblue')
rects2 = ax.bar(x + width/2, avg_similarities_str, width, yerr=std_dev_str, label='Average similarity to Str', capsize=2, color='salmon')

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Average Cosine Similarity')
ax.set_xticks(x)
ax.set_xticklabels(str_labels)
ax.legend()

plt.xticks(rotation=45, ha='right')
plt.tight_layout()

# Save the plot as an SVG file
plt.savefig('average_cosine_similarity_with_error_bars.svg')

# Return the figure for display
plt.show()