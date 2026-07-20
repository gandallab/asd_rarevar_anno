import scanpy as sc
from anndata import read_h5ad
import pandas as pd
import numpy as np
import os
from os.path import join
import time
import argparse
import matplotlib.pyplot as plt
from scipy import sparse

#  autoreload
%load_ext autoreload
%autoreload 2

# scMultiome 
H5AD_FILE='/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/obj_rna_raw.h5ad'
adata = read_h5ad(H5AD_FILE)

# cov 
df_cov = pd.DataFrame(index=adata.obs.index)
df_cov['const'] = 1 # constant of cov
df_cov['n_genes'] = adata.obs['nFeature_RNA']
df_cov['sex_male'] = (adata.obs['sex']=='M')*1

# create dummy cov variable for donor ID
for donor in sorted(set(adata.obs['Ident'])):
        df_cov['donor_%s'%donor] = (adata.obs['Ident']==donor)*1

df_cov.to_csv('/mnt/isilon/gandal_lab/liaoyd/project/asd_rarevar_anno/data/SnMultiome_Wang2025/donorID_sex_ngene.cov',sep='\t')
