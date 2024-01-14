# peptide_similarity
3 methods to find peptide similarity:

cosine.py gives the cosine similarity and the cosine_similarities.txt and average_cosine_similarity_with_error_bars.svg files.

map4_similarity.py gives the Jaccard distance (with MAP4 fingerprints) using the whole peptide SMILES. It uses the smiles.smi file (with the Pep1-12 and Str1-9 written in SMILES) and gives the minhash_distance.txt and minhash_distance_similarity.svg files.

map4_similarity_positional.py gives the Jaccard distance (with MAP4 fingerprints) using each amino acid individually and by position. It uses the converted_sequences.csv file (with a list of the amino acids individually written for each position for Pep1-12 and Str1-9) and gives the amino_acid_position_similarity.csv and minhash_distance_similarity_12aminoacids.svg files.
