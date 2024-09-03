# monae
[![DOI](https://zenodo.org/badge/830540235.svg)](https://zenodo.org/doi/10.5281/zenodo.13636951)

monae: multi-modal single-cell integration and imputation

# Environment
Python version: python 3.8

Install the following dependencies via pip or conda:
```python
pip install muon
pip install git+https://github.com/pinellolab/simba_pbg.git
pip install git+https://github.com/huidongchen/simba
pip install scanpy
pip install scglue
pip install pytorch-ignite
pip install typing_extensions
pip install tensorboardX
pip install torch-1.11.0+cu113-cp38-cp38-linux_x86_64.whl
```

# Dataset
Dataset processed (`dataset_monae`)：

link: https://pan.baidu.com/s/1lOjSVVIUQ3WmqTDUJzdxMQ?pwd=1etw 

For the raw data, select annotation file in `preprocess_1.py`:
```python
gtf="./src/bedtools/annotation/gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz", # Human
gtf="./src/bedtools/annotation/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz", # Mouse
```

# Example
See the monae example in the `notebook_monae`. Trained Monae model:

link: https://pan.baidu.com/s/1F-EzYazz31_LEOMEfxMM2w?pwd=v8wl

The models and datasets are organized as follows:
```python
--src
...
--notebook_monae
  --dataset_monae
  --model # trained monae model
  --integration.ipynb # for integration
  --imputation_rna2atac.ipynb # for cross-modal imputation (rna->atac)
  --imputation_atac2rna.ipynb # for cross-modal imputation (atac->rna)
```

In ipynb files, select one of four datasets: `pbmc, chen2019, ma2020, muto2021`.

# Monae-Extension
A more stable and faster training strategy is provided here for Monae's variant Monae-E. Run the following files in sequence:

1. preprocess_1.py: data preprocess and constructing graph guidance.
2. preprocess_2.py: Unified learning of cell and feature representations.
3. integration.py: Train model.

Take the PBMC dataset as an example, the trained model has been saved here: `dataset/pbmc/result/ckpt/monae.dill`. Loading trained model in `inference.ipynb`.

# reference

[1] Cellcano: supervised cell type identification for single cell ATAC-seq data, nature communications, 2023

https://github.com/marvinquiet/Cellcano

[2] Multi-omics single-cell data integration and regulatory inference with graph-linked embedding, nature biotechnology, 2022

https://github.com/gao-lab/GLUE

[3] SIMBA: single-cell embedding along with features, nature methods, 2023

https://github.com/pinellolab/simba

[4] scGPT: toward building a foundation model for single-cell multi-omics using generative AI, nature methods, 2024

https://github.com/bowang-lab/scGPT
