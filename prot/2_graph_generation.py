import os
import sys
sys.path.append("..")
import src.emb_init as si
import anndata
import matplotlib.pyplot as plt
import time 
import muon as mu
import scanpy as sc
from scipy.sparse import csr_matrix
import networkx as nx
import itertools
import scglue
import warnings
warnings.filterwarnings("ignore")

workdir = f"./gnn_result"

mdata = mu.read_h5mu("./adatas/opcite_preprocessed.h5mu.gz")
print(mdata)
adata_CPro = mdata['adt']
adata_CG = mdata['rna']

# ATAC预处理
# si.pp.cal_qc_rna(adata_CPro)
si.pp.normalize(adata_CPro,method='lib_size')
si.pp.log_transform(adata_CPro)

# RNA预处理
# si.pp.cal_qc_rna(adata_CG)
si.pp.normalize(adata_CG,method='lib_size')
si.pp.log_transform(adata_CG)

# print(adata_CG.var_names)
# print(adata_CP.var_names)

# load in graph ('graph0') info
si.load_graph_stats(path=f'{workdir}/pbg/graph0/')
# load in model info for ('graph0')
si.load_pbg_config(path=f'{workdir}/pbg/graph0/model/')

dict_adata = si.read_embedding()
adata_C = dict_adata['C']  # embeddings for cells
adata_G = dict_adata['G']  # embeddings for genes
adata_Pro = dict_adata['Pro']  # embeddings for peaks

# print(adata_G.obs_names)
# print(adata_P.obs_names)

adata_all = si.tl.embed(adata_ref=adata_C, list_adata_query=[adata_G, adata_Pro])
# adata_CP.varm["X_node"] = adata_all[adata_CP.var_names,:].X
# adata_CG.varm["X_node"] = adata_all[adata_CG.var_names,:].X

adata_G.X = adata_all[adata_G.obs_names,:].X
adata_Pro.X = adata_all[adata_Pro.obs_names,:].X

adata_CrnaCatac = si.tl.infer_guidance_edges(adata_G, adata_Pro, n_components=10, k=36)
# print(type(adata_CrnaCatac.layers['conn']))
# print(adata_CrnaCatac.layers['conn'].shape)
# print(adata_CrnaCatac.obs_names[0])
# print(adata_CrnaCatac.var_names[1])


graph = nx.MultiDiGraph()
for i in range(adata_CrnaCatac.layers['conn'].shape[0]):
    for j in range(adata_CrnaCatac.layers['conn'].shape[1]):
        if adata_CrnaCatac.layers['conn'][i, j] != 0:
            graph.add_edge(adata_CrnaCatac.obs_names[i], adata_CrnaCatac.var_names[j], weight=1.0, sign=1)

graph = scglue.graph.compose_multigraph(graph, graph.reverse())
for i in itertools.chain(adata_CG.obs_names, adata_CPro.var_names):
    graph.add_edge(i, i, weight=1.0, sign=1)

print(graph)
print("Nodes:")
for node in list(graph.nodes())[:10]:
    print(f"Node: {node}, Attributes: {graph.nodes[node]}")
print("\nEdges:")
for edge in list(graph.edges(data=True))[:10]:
    print(f"Edge: {edge[0]} -> {edge[1]}, Attributes: {edge[2]}")

nx.write_graphml(graph, "./adatas/guidance-hvf.graphml.gz")

adata_CG.obsm['X_init'] = adata_all[adata_CG.obs_names,:].X
adata_CPro.obsm['X_init'] = adata_all[adata_CPro.obs_names,:].X

print(adata_CG.X.max())
print(adata_CPro.X.max())

adata_CG.write("./adatas/rna_hvg.h5ad", compression="gzip")
adata_CPro.write("./adatas/adt_hvg.h5ad", compression="gzip")
