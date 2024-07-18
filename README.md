# monae
monae: multi-modal single-cell integration and imputation

# Environment
Python version: python 3.8

Install the following dependencies via pip or conda:
```python
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
Dataset processed：

link：https://pan.baidu.com/s/15Jb_KwllAUaGAOm11bJFbg?pwd=tw6z

For the raw data, select annotation file in `preprocess_1.py`:
```python
gtf="./src/bedtools/annotation/gencode.v42.chr_patch_hapl_scaff.annotation.gtf.gz", # Human
gtf="./src/bedtools/annotation/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf.gz", # Mouse
```


# reference

[1] Cellcano: supervised cell type identification for single cell ATAC-seq data, nature communications, 2023

https://github.com/marvinquiet/Cellcano

[2] Multi-omics single-cell data integration and regulatory inference with graph-linked embedding, nature biotechnology, 2022

https://github.com/gao-lab/GLUE

[3] SIMBA: single-cell embedding along with features, nature methods, 2023

https://github.com/pinellolab/simba

[4] scGPT: toward building a foundation model for single-cell multi-omics using generative AI, nature methods, 2024

https://github.com/bowang-lab/scGPT
