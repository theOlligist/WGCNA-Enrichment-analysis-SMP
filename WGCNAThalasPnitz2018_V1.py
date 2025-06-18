# metatranscriptomics_pipeline.py
# Python translation of WGCNA + edgeR R pipeline (approximate)

import pandas as pd
import numpy as np
import re
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster
from scipy.spatial.distance import pdist
from scipy.stats import pearsonr
import seaborn as sns
import matplotlib.pyplot as plt

# -------------------- 1. Load & preprocess raw count table --------------------
# Simulate Diatom_raw_df equivalent
# Placeholder: replace with real file loading
Diatom_raw_df = pd.DataFrame()  # load data here

# Filter by taxonomy and add custom label
Diatom_raw_df = Diatom_raw_df[Diatom_raw_df['Taxonomy'].str.contains("Thalassiosirales|Bacillariales")].copy()
Diatom_raw_df['tax'] = Diatom_raw_df['Taxonomy'].apply(
    lambda x: 'Thalassiosira' if 'Thalassiosirales' in x else 'Pseudo_nitz')
Diatom_raw_df['ID'] = Diatom_raw_df['tax'] + '-' + Diatom_raw_df['KO']

# Aggregate counts by ID and sample
value_vars = [col for col in Diatom_raw_df.columns if col.startswith("SMP")]
long_df = Diatom_raw_df.melt(id_vars=['ID'], value_vars=value_vars, var_name='sample', value_name='count')
long_df['sample'] = long_df['sample'].str.replace(r"SMPier\\.|\\.\\d{2}\\.\\d{2}\\.2018_S\\d+", "", regex=True)

agg_df = long_df.groupby(['ID', 'sample'])['count'].sum().reset_index()
wide_df = agg_df.pivot(index='ID', columns='sample', values='count').fillna(0)

# -------------------- 2. Normalize with CPM (approx. TMM) --------------------
cpm_df = wide_df.div(wide_df.sum(axis=0), axis=1) * 1e6

# Round and log transform
Norm_df = np.round(np.log2(cpm_df + 1))  # log2(count + 1)
Norm_df_t = Norm_df.transpose()

# -------------------- 3. Soft-threshold and Adjacency matrix --------------------
def adjacency_matrix(data, power=6):
    corr_matrix = data.corr(method='pearson')
    return np.power(np.abs(corr_matrix), power)

powers = list(range(1, 11)) + list(range(12, 21, 2))

# Placeholder: you would evaluate soft threshold fit here
soft_power = 6
adj_mat = adjacency_matrix(Norm_df_t, power=soft_power)

# -------------------- 4. TOM (Topological Overlap) --------------------
def TOM_similarity(adj):
    n = adj.shape[0]
    TOM = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i == j:
                TOM[i, j] = 1
            else:
                numer = np.sum(np.minimum(adj.iloc[i, :], adj.iloc[j, :])) + adj.iloc[i, j]
                denom = min(np.sum(adj.iloc[i, :]), np.sum(adj.iloc[j, :])) + 1 - adj.iloc[i, j]
                TOM[i, j] = numer / denom
    return pd.DataFrame(TOM, index=adj.index, columns=adj.columns)

TOM_df = TOM_similarity(adj_mat)
dissTOM = 1 - TOM_df

# -------------------- 5. Hierarchical clustering & module detection --------------------
dist_matrix = pdist(dissTOM)
linkage_matrix = linkage(dist_matrix, method='average')

# Mimic dynamic tree cut by choosing distance threshold manually
module_labels = fcluster(linkage_matrix, t=0.25, criterion='distance')

# Assign colors (optional)
module_colors = pd.Series(module_labels, index=Norm_df.index).map(lambda x: f"module_{x}")

# -------------------- 6. Module Eigengenes --------------------
eigengenes = {}
for label in module_colors.unique():
    genes = module_colors[module_colors == label].index
    pca = PCA(n_components=1)
    eigengene = pca.fit_transform(Norm_df.loc[genes].T)
    eigengenes[label] = eigengene.flatten()
ME_df = pd.DataFrame(eigengenes, index=Norm_df.columns)

# -------------------- 7. Trait correlation --------------------
# Placeholder: replace with real environmental data
trait_df = pd.DataFrame()  # rows = samples, cols = traits
cor_matrix = ME_df.corrwith(trait_df, axis=0)

# -------------------- 8. Enrichment analysis --------------------
# Optional with gseapy or bioservices; placeholders below
# import gseapy as gp
# enr = gp.enrichr(gene_list=[...], gene_sets='KEGG_2019_Human', organism='Human')

# -------------------- 9. DE analysis (approximate) --------------------
# Use statsmodels or scipy to compare groups manually; edgeR not replicated directly

# -------------------- 10. Plotting (optional) --------------------
# Example dendrogram:
# dendrogram(linkage_matrix, labels=Norm_df.index)

print("Pipeline executed (placeholder mode). Replace mock inputs with real data.")
