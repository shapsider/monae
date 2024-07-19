import anndata
import networkx as nx
import scanpy as sc
from anndata import AnnData
import torch
from typing import Any, List, Mapping, Optional
import pandas as pd

class SingletonMeta(type):

    r"""
    Ensure singletons via a meta class
    """

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super().__call__(*args, **kwargs)
        return cls._instances[cls]

class ConfigManager(metaclass=SingletonMeta):

    r"""
    Global configurations
    """

    def __init__(self) -> None:
        self.TMP_PREFIX = "COVELTMP"
        self.ANNDATA_KEY = "__covel__"
        self.CPU_ONLY = False
        self.CUDNN_MODE = "repeatability"
        self.MASKED_GPUS = []
        self.ARRAY_SHUFFLE_NUM_WORKERS = 0
        self.GRAPH_SHUFFLE_NUM_WORKERS = 1
        self.FORCE_TERMINATE_WORKER_PATIENCE = 60
        self.DATALOADER_NUM_WORKERS = 0
        self.DATALOADER_FETCHES_PER_WORKER = 4
        self.DATALOADER_PIN_MEMORY = True
        self.CHECKPOINT_SAVE_INTERVAL = 10
        self.CHECKPOINT_SAVE_NUMBERS = 3
        self.PRINT_LOSS_INTERVAL = 10
        self.TENSORBOARD_FLUSH_SECS = 5
        self.ALLOW_TRAINING_INTERRUPTION = True

    @property
    def TMP_PREFIX(self) -> str:
        r"""
        Prefix of temporary files and directories created.
        Default values is ``"COVELTMP"``.
        """
        return self._TMP_PREFIX

    @TMP_PREFIX.setter
    def TMP_PREFIX(self, tmp_prefix: str) -> None:
        self._TMP_PREFIX = tmp_prefix

    @property
    def ANNDATA_KEY(self) -> str:
        r"""
        Key in ``adata.uns`` for storing dataset configurations.
        Default value is ``"__covel__"``
        """
        return self._ANNDATA_KEY

    @ANNDATA_KEY.setter
    def ANNDATA_KEY(self, anndata_key: str) -> None:
        self._ANNDATA_KEY = anndata_key

    @property
    def CPU_ONLY(self) -> bool:
        r"""
        Whether computation should use only CPUs.
        Default value is ``False``.
        """
        return self._CPU_ONLY

    @CPU_ONLY.setter
    def CPU_ONLY(self, cpu_only: bool) -> None:
        self._CPU_ONLY = cpu_only
        if self._CPU_ONLY and self._DATALOADER_NUM_WORKERS:
            self.logger.warning(
                "It is recommended to set `DATALOADER_NUM_WORKERS` to 0 "
                "when using CPU_ONLY mode. Otherwise, deadlocks may happen "
                "occationally."
            )

    @property
    def CUDNN_MODE(self) -> str:
        r"""
        CuDNN computation mode, should be one of {"repeatability", "performance"}.
        Default value is ``"repeatability"``.

        Note
        ----
        As of now, due to the use of :meth:`torch.Tensor.scatter_add_`
        operation, the results are not completely reproducible even when
        ``CUDNN_MODE`` is set to ``"repeatability"``, if GPU is used as
        computation device. Exact repeatability can only be achieved on CPU.
        The situtation might change with new releases of :mod:`torch`.
        """
        return self._CUDNN_MODE

    @CUDNN_MODE.setter
    def CUDNN_MODE(self, cudnn_mode: str) -> None:
        if cudnn_mode not in ("repeatability", "performance"):
            raise ValueError("Invalid mode!")
        self._CUDNN_MODE = cudnn_mode
        torch.backends.cudnn.deterministic = self._CUDNN_MODE == "repeatability"
        torch.backends.cudnn.benchmark = self._CUDNN_MODE == "performance"

    @property
    def MASKED_GPUS(self) -> List[int]:
        r"""
        A list of GPUs that should not be used when selecting computation device.
        This must be set before initializing any model, otherwise would be ineffective.
        Default value is ``[]``.
        """
        return self._MASKED_GPUS

    @MASKED_GPUS.setter
    def MASKED_GPUS(self, masked_gpus: List[int]) -> None:
        if masked_gpus:
            import pynvml
            pynvml.nvmlInit()
            device_count = pynvml.nvmlDeviceGetCount()
            for item in masked_gpus:
                if item >= device_count:
                    raise ValueError(f"GPU device \"{item}\" is non-existent!")
        self._MASKED_GPUS = masked_gpus

    @property
    def ARRAY_SHUFFLE_NUM_WORKERS(self) -> int:
        r"""
        Number of background workers for array data shuffling.
        Default value is ``0``.
        """
        return self._ARRAY_SHUFFLE_NUM_WORKERS

    @ARRAY_SHUFFLE_NUM_WORKERS.setter
    def ARRAY_SHUFFLE_NUM_WORKERS(self, array_shuffle_num_workers: int) -> None:
        self._ARRAY_SHUFFLE_NUM_WORKERS = array_shuffle_num_workers

    @property
    def GRAPH_SHUFFLE_NUM_WORKERS(self) -> int:
        r"""
        Number of background workers for graph data shuffling.
        Default value is ``1``.
        """
        return self._GRAPH_SHUFFLE_NUM_WORKERS

    @GRAPH_SHUFFLE_NUM_WORKERS.setter
    def GRAPH_SHUFFLE_NUM_WORKERS(self, graph_shuffle_num_workers: int) -> None:
        self._GRAPH_SHUFFLE_NUM_WORKERS = graph_shuffle_num_workers

    @property
    def FORCE_TERMINATE_WORKER_PATIENCE(self) -> int:
        r"""
        Seconds to wait before force terminating unresponsive workers.
        Default value is ``60``.
        """
        return self._FORCE_TERMINATE_WORKER_PATIENCE

    @FORCE_TERMINATE_WORKER_PATIENCE.setter
    def FORCE_TERMINATE_WORKER_PATIENCE(self, force_terminate_worker_patience: int) -> None:
        self._FORCE_TERMINATE_WORKER_PATIENCE = force_terminate_worker_patience

    @property
    def DATALOADER_NUM_WORKERS(self) -> int:
        r"""
        Number of worker processes to use in data loader.
        Default value is ``0``.
        """
        return self._DATALOADER_NUM_WORKERS

    @DATALOADER_NUM_WORKERS.setter
    def DATALOADER_NUM_WORKERS(self, dataloader_num_workers: int) -> None:
        if dataloader_num_workers > 8:
            self.logger.warning(
                "Worker number 1-8 is generally sufficient, "
                "too many workers might have negative impact on speed."
            )
        self._DATALOADER_NUM_WORKERS = dataloader_num_workers

    @property
    def DATALOADER_FETCHES_PER_WORKER(self) -> int:
        r"""
        Number of fetches per worker per batch to use in data loader.
        Default value is ``4``.
        """
        return self._DATALOADER_FETCHES_PER_WORKER

    @DATALOADER_FETCHES_PER_WORKER.setter
    def DATALOADER_FETCHES_PER_WORKER(self, dataloader_fetches_per_worker: int) -> None:
        self._DATALOADER_FETCHES_PER_WORKER = dataloader_fetches_per_worker

    @property
    def DATALOADER_FETCHES_PER_BATCH(self) -> int:
        r"""
        Number of fetches per batch in data loader (read-only).
        """
        return max(1, self.DATALOADER_NUM_WORKERS) * self.DATALOADER_FETCHES_PER_WORKER

    @property
    def DATALOADER_PIN_MEMORY(self) -> bool:
        r"""
        Whether to use pin memory in data loader.
        Default value is ``True``.
        """
        return self._DATALOADER_PIN_MEMORY

    @DATALOADER_PIN_MEMORY.setter
    def DATALOADER_PIN_MEMORY(self, dataloader_pin_memory: bool):
        self._DATALOADER_PIN_MEMORY = dataloader_pin_memory

    @property
    def CHECKPOINT_SAVE_INTERVAL(self) -> int:
        r"""
        Automatically save checkpoints every n epochs.
        Default value is ``10``.
        """
        return self._CHECKPOINT_SAVE_INTERVAL

    @CHECKPOINT_SAVE_INTERVAL.setter
    def CHECKPOINT_SAVE_INTERVAL(self, checkpoint_save_interval: int) -> None:
        self._CHECKPOINT_SAVE_INTERVAL = checkpoint_save_interval

    @property
    def CHECKPOINT_SAVE_NUMBERS(self) -> int:
        r"""
        Maximal number of checkpoints to preserve at any point.
        Default value is ``3``.
        """
        return self._CHECKPOINT_SAVE_NUMBERS

    @CHECKPOINT_SAVE_NUMBERS.setter
    def CHECKPOINT_SAVE_NUMBERS(self, checkpoint_save_numbers: int) -> None:
        self._CHECKPOINT_SAVE_NUMBERS = checkpoint_save_numbers

    @property
    def PRINT_LOSS_INTERVAL(self) -> int:
        r"""
        Print loss values every n epochs.
        Default value is ``10``.
        """
        return self._PRINT_LOSS_INTERVAL

    @PRINT_LOSS_INTERVAL.setter
    def PRINT_LOSS_INTERVAL(self, print_loss_interval: int) -> None:
        self._PRINT_LOSS_INTERVAL = print_loss_interval

    @property
    def TENSORBOARD_FLUSH_SECS(self) -> int:
        r"""
        Flush tensorboard logs to file every n seconds.
        Default values is ``5``.
        """
        return self._TENSORBOARD_FLUSH_SECS

    @TENSORBOARD_FLUSH_SECS.setter
    def TENSORBOARD_FLUSH_SECS(self, tensorboard_flush_secs: int) -> None:
        self._TENSORBOARD_FLUSH_SECS = tensorboard_flush_secs

    @property
    def ALLOW_TRAINING_INTERRUPTION(self) -> bool:
        r"""
        Allow interruption before model training converges.
        Default values is ``True``.
        """
        return self._ALLOW_TRAINING_INTERRUPTION

    @ALLOW_TRAINING_INTERRUPTION.setter
    def ALLOW_TRAINING_INTERRUPTION(self, allow_training_interruption: bool) -> None:
        self._ALLOW_TRAINING_INTERRUPTION = allow_training_interruption


config = ConfigManager()

def configure_dataset(
        adata: AnnData, prob_model: str,
        use_highly_variable: bool = True,
        use_layer: Optional[str] = None,
        use_rep: Optional[str] = None,
        use_batch: Optional[str] = None,
        use_cell_type: Optional[str] = None,
        use_dsc_weight: Optional[str] = None,
        use_obs_names: bool = False
) -> None:
    data_config = {}
    data_config["prob_model"] = prob_model
    if use_highly_variable:
        if "highly_variable" not in adata.var:
            raise ValueError("Please mark highly variable features first!")
        data_config["use_highly_variable"] = True
        data_config["features"] = adata.var.query("highly_variable").index.to_numpy().tolist()
    else:
        data_config["use_highly_variable"] = False
        data_config["features"] = adata.var_names.to_numpy().tolist()
    if use_layer:
        if use_layer not in adata.layers:
            raise ValueError("Invalid `use_layer`!")
        data_config["use_layer"] = use_layer
    else:
        data_config["use_layer"] = None
    if use_rep:
        if use_rep not in adata.obsm:
            raise ValueError("Invalid `use_rep`!")
        data_config["use_rep"] = use_rep
        data_config["rep_dim"] = adata.obsm[use_rep].shape[1]
    else:
        data_config["use_rep"] = None
        data_config["rep_dim"] = None
    if use_batch:
        if use_batch not in adata.obs:
            raise ValueError("Invalid `use_batch`!")
        data_config["use_batch"] = use_batch
        data_config["batches"] = pd.Index(
            adata.obs[use_batch]
        ).dropna().drop_duplicates().sort_values().to_numpy()  # AnnData does not support saving pd.Index in uns
    else:
        data_config["use_batch"] = None
        data_config["batches"] = None
    if use_cell_type:
        if use_cell_type not in adata.obs:
            raise ValueError("Invalid `use_cell_type`!")
        data_config["use_cell_type"] = use_cell_type
        data_config["cell_types"] = pd.Index(
            adata.obs[use_cell_type]
        ).dropna().drop_duplicates().sort_values().to_numpy()  # AnnData does not support saving pd.Index in uns
    else:
        data_config["use_cell_type"] = None
        data_config["cell_types"] = None
    if use_dsc_weight:
        if use_dsc_weight not in adata.obs:
            raise ValueError("Invalid `use_dsc_weight`!")
        data_config["use_dsc_weight"] = use_dsc_weight
    else:
        data_config["use_dsc_weight"] = None
    data_config["use_obs_names"] = use_obs_names
    adata.uns[config.ANNDATA_KEY] = data_config


if __name__=="__main__":
    print(config.ANNDATA_KEY)
    atac = anndata.read_h5ad("../data/atac_hvg.h5ad")
    configure_dataset(atac, "NB", use_highly_variable=True, use_rep="pre")
    print(atac.uns[config.ANNDATA_KEY])