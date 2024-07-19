import anndata
import networkx as nx
import scanpy as sc
from anndata import AnnData
import torch
from typing import Any, List, Mapping, Optional
import pandas as pd
import os
import numpy as np

from pathlib import Path

from .typehint import Kws
from .models import CoVELModel
from .config import config
from .auxtools import get_metacells
from .config import configure_dataset

def covel_train(
        adatas: List[AnnData], graph: nx.Graph, 
        model: type = CoVELModel,
        init_kws: Kws = None, compile_kws: Kws = None, fit_kws: Kws = None,
        balance_kws: Kws = None,
        config = None,
        result_path = None
) -> CoVELModel:
    [modal_names, prob, rep] = config

    for idx, adata in enumerate(adatas):
        configure_dataset(adata, prob[idx], 
                      use_highly_variable=True,
                      use_rep=rep[idx])

    adatas = dict(zip(modal_names, adatas))

    init_kws = init_kws or {}
    compile_kws = compile_kws or {}
    fit_kws = fit_kws or {}
    balance_kws = balance_kws or {}

    print("Training...")
    finetune_fit_kws = fit_kws.copy()

    finetune = CoVELModel(adatas, sorted(graph.nodes), **init_kws)
    finetune.compile(**compile_kws)
    finetune.fit(adatas, graph, **finetune_fit_kws)
    if "directory" in finetune_fit_kws:
        finetune.save(os.path.join(finetune_fit_kws["directory"], "monae.dill"))
    
    return finetune