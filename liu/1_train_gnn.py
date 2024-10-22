import os
import sys
sys.path.append("..")
import src.emb_init as si
import anndata
import matplotlib.pyplot as plt
import time 
import muon as mu
import scanpy as sc
import numpy as np
import warnings
import torch
warnings.filterwarnings("ignore")

def find_indexes(list1, list2):
    indexes = []
    for item in list2:
        index = list1.index(item)
        indexes.append(index)
    return indexes

def generate_gaussian_vectors(lst, length, mean=0, std_dev=1):
    gaussian_vectors = []
    for _ in lst:
        vector = np.random.normal(mean, std_dev, size=(length,))
        gaussian_vectors.append(vector)
    return np.array(gaussian_vectors)

# Load data into a Muon object.
mdata = mu.read_h5mu("./adatas/liu_preprocessed.h5mu.gz")
print(mdata)

adata_CP = mdata['atac']
adata_CG = mdata['rna']

print(adata_CG.X.max())
print(adata_CP.X.max())

workdir = f"./gnn_result"
si.settings.set_workdir(workdir)

# ATAC预处理
# si.pp.filter_genes(adata_CP,min_n_cells=3)
si.pp.cal_qc_rna(adata_CP)
si.pp.normalize(adata_CP,method='lib_size')
si.pp.log_transform(adata_CP)

si.tl.discretize(adata_CP,n_bins=5)
si.pl.discretize(adata_CP,kde=False)
plt.savefig(f"{workdir}/discretize_atac.png")

# RNA预处理
# si.pp.filter_genes(adata_CG,min_n_cells=3)
si.pp.cal_qc_rna(adata_CG)
si.pp.normalize(adata_CG,method='lib_size')
si.pp.log_transform(adata_CG)

si.tl.discretize(adata_CG,n_bins=5)
si.pl.discretize(adata_CG,kde=False)
plt.savefig(f"{workdir}/discretize_rna.png")

# 生成图
start_time = time.time()
si.tl.gen_graph(list_CP=[adata_CP],
                list_CG=[adata_CG],
                copy=False,
                use_highly_variable=False,
                use_top_pcs=False,
                dirname='graph0')
end_time = time.time()
run_time = end_time - start_time
print(f"run time: {run_time}s")

if __name__ == '__main__':
    start_time = time.time()
    si.tl.pbg_train(auto_wd=True, 
                    save_wd=True, 
                    output='model', 
                    pre_embeddings=None,
                    )
    end_time = time.time()
    run_time = end_time - start_time
    print(f"run time: {run_time}s")
    si.pl.pbg_metrics(fig_ncol=3)
    plt.savefig(f"{workdir}/pbg_metrics.png")
