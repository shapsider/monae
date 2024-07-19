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

from src.config import configure_dataset
from src.train import covel_train

dataset = "pbmc"
result_path = f"./dataset/{dataset}/result"
if(not os.path.exists(result_path)):
            os.makedirs(result_path)

rna = anndata.read_h5ad(f"./dataset/{dataset}/rna_hvg.h5ad")
atac = anndata.read_h5ad(f"./dataset/{dataset}/atac_hvg.h5ad")
graph = nx.read_graphml(f"./dataset/{dataset}/guidance-hvf.graphml.gz")

rna.obs['cell type'] = rna.obs['cell_type']
atac.obs['cell type'] = atac.obs['cell_type']

rna.X=rna.layers['counts']
atac.X=atac.layers['counts']

adatas=[rna, atac]
modal_names=["RNA", "ATAC"]
prob=['NB','NB']
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