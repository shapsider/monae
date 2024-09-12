import anndata
import networkx as nx
import scanpy as sc
from anndata import AnnData
import pandas as pd
import os
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import time
import sys
sys.path.append("..")
from src.config import configure_dataset
from src.train import covel_train


result_path = f"./result"
if(not os.path.exists(result_path)):
            os.makedirs(result_path)

rna = anndata.read_h5ad(f"./adatas/rna_hvg.h5ad")
atac = anndata.read_h5ad(f"./adatas/atac_hvg.h5ad")
graph = nx.read_graphml(f"./adatas/guidance-hvf.graphml.gz")

rna.obs['cell type'] = rna.obs['celltype']
atac.obs['cell type'] = atac.obs['celltype']

adatas=[rna, atac]
modal_names=["RNA", "ATAC"]
prob=['Normal','Normal']
rep = ['X_init', 'X_init']
save_path = f"{result_path}/ckpt"

vertices = sorted(graph.nodes)

for idx, adata in enumerate(adatas):
    configure_dataset(adata, prob[idx], 
                      use_highly_variable=True,
                      use_rep=rep[idx],
                      )

data = dict(zip(modal_names, adatas))

starttime = time.time()
covel = covel_train(
        adatas, 
        graph,
        fit_kws={"directory": save_path},
        config = [modal_names, prob, rep],
        result_path = result_path
)
endtime = time.time()
print("#########################################################################")
print('spent {:.5f} s'.format(endtime-starttime))
print("#########################################################################")