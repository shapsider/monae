import anndata as ad
import networkx as nx
import scanpy as sc
import scglue
import numpy as np
from src.bedtools import genomics

import matplotlib.pyplot as plt

dataset = "pbmc"
rna = ad.read_h5ad(f"./dataset/{dataset}/raw/10x-Multiome-Pbmc10k-RNA.h5ad")
atac = ad.read_h5ad(f"./dataset/{dataset}/raw/10x-Multiome-Pbmc10k-ATAC.h5ad")
print(rna)
print(atac)

sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor="seurat_v3")
rna_hvg = rna[:, rna.var.highly_variable]
rna_hvg.layers["counts"] = rna_hvg.X.copy()
sc.pp.normalize_total(rna_hvg)
sc.pp.log1p(rna_hvg)
sc.pp.scale(rna_hvg)

sc.tl.pca(rna_hvg, n_comps=10, svd_solver="auto")
sc.pp.neighbors(rna_hvg, use_rep="X_pca", metric="cosine")
sc.tl.umap(rna_hvg)
sc.pl.umap(rna_hvg, color="cell_type", legend_loc='on data')

# TODO
scglue.data.get_gene_annotation(
    rna, 
    gtf="./src/bedtools/annotation/gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz",
    # gtf="./src/bedtools/annotation/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz",
    gtf_by="gene_name"
)

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)

# graph
guidance = genomics.rna_anchored_guidance_graph(rna, atac)
print(guidance)
scglue.graph.check_graph(guidance, [rna, atac])
rna_hvg.write(f"./dataset/{dataset}/rna_hvg.h5ad", compression="gzip")

# atac-hvg
hvg_reachable = scglue.graph.reachable_vertices(guidance, rna.var.query("highly_variable").index)

atac.var["highly_variable"] = [item in hvg_reachable for item in atac.var_names]
print(atac.var["highly_variable"].sum())

atac_hvg = atac[:, atac.var.highly_variable]
atac_hvg.layers["counts"] = atac_hvg.X.copy()
sc.pp.normalize_total(atac_hvg)
sc.pp.log1p(atac_hvg)
sc.pp.scale(atac_hvg)
print(atac_hvg)

sc.tl.pca(atac_hvg, n_comps=10, svd_solver="auto")
sc.pp.neighbors(atac_hvg, use_rep="X_pca", metric="cosine")
sc.tl.umap(atac_hvg)
sc.pl.umap(atac_hvg, color="cell_type")

atac_hvg.write(f"./dataset/{dataset}/atac_hvg.h5ad", compression="gzip")

subgraph = guidance.subgraph(hvg_reachable)
nx.write_graphml(subgraph, f"./dataset/{dataset}/guidance-hvf.graphml.gz")

print(rna_hvg.layers['counts'].max())
print(atac_hvg.layers['counts'].max())