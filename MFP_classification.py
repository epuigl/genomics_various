# Python3.7

import pandas as pd
import numpy as np
import importlib.util

#---- Import custom modules
spec = importlib.util.spec_from_file_location("MFP","MFP/portraits/utils.py")
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)

spec = importlib.util.spec_from_file_location("MFP","MFP/portraits/clustering.py")
clustering = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clustering)

spec = importlib.util.spec_from_file_location("MFP","MFP/portraits/classification.py")
classification = importlib.util.module_from_spec(spec)
spec.loader.exec_module(classification)

#---- Classify another cohort having clusters on TCGA as training

# Load training cohort with known MFP labels
TCGA_signature_scores_scaled = pd.read_csv("MFP/Cohorts/Pan_TCGA/signatures.tsv", sep='\t', index_col=0).T #signatures in rows

TCGA_annotation = pd.read_csv("MFP/Cohorts/Pan_TCGA/annotation.tsv", sep='\t', index_col=0)

# Fit the model
MODEL = classification.KNeighborsClusterClassifier(norm=False, scale=False, clip=2, k=35).fit(TCGA_signature_scores_scaled,
                                                   TCGA_annotation.MFP)

# Load the cohort of interest
# Read signatures
gmt = utils.read_gene_sets('MFP/signatures/gene_signatures.gmt') #GMT format like in MSIGdb

# Read expressions
exp = pd.read_csv("my_counts.tsv", sep='\t', header=0, index_col=0).T  #TPM/CPM values; genes in rows
    
if exp.max().max() > 35:
    print('Performing log2+1 transformation')
    exp = np.log2(1+exp)

# Calculate signature scores for the cohort
signature_scores = utils.ssgsea_formula(exp, gmt)

# Scale signatures
signature_scores_scaled = utils.median_scale(signature_scores, 2)

# Predict clusters
cluster_labels = MODEL.predict(signature_scores_scaled[MODEL.X.columns]).rename('MFP')
print('Predicted labels count:')
print(cluster_labels.value_counts())

# Output the clusters
cluster_labels.to_csv("MFP_labels_cohort.tsv", sep='\t', index=True)