r"""
Miscellaneous utilities
"""

import os
import logging
import signal
import subprocess
import sys
from collections import defaultdict
from multiprocessing import Process
from typing import Any, List, Mapping, Optional

import numpy as np
import pandas as pd

from .typehint import RandomState, T

from .config import SingletonMeta

AUTO = "AUTO"  # Flag for using automatically determined hyperparameters


#------------------------------ Global containers ------------------------------

processes: Mapping[int, Mapping[int, Process]] = defaultdict(dict)  # id -> pid -> process

#--------------------------------- Log manager ---------------------------------

class _CriticalFilter(logging.Filter):

    def filter(self, record: logging.LogRecord) -> bool:
        return record.levelno >= logging.WARNING


class _NonCriticalFilter(logging.Filter):

    def filter(self, record: logging.LogRecord) -> bool:
        return record.levelno < logging.WARNING


class LogManager(metaclass=SingletonMeta):

    r"""
    Manage loggers used in the package
    """

    def __init__(self) -> None:
        self._loggers = {}
        self._log_file = None
        self._console_log_level = logging.INFO
        self._file_log_level = logging.DEBUG
        self._file_fmt = \
            "%(asctime)s.%(msecs)03d [%(levelname)s] %(name)s: %(message)s"
        self._console_fmt = \
            "[%(levelname)s] %(name)s: %(message)s"
        self._date_fmt = "%Y-%m-%d %H:%M:%S"

    @property
    def log_file(self) -> str:
        r"""
        Configure log file
        """
        return self._log_file

    @property
    def file_log_level(self) -> int:
        r"""
        Configure logging level in the log file
        """
        return self._file_log_level

    @property
    def console_log_level(self) -> int:
        r"""
        Configure logging level printed in the console
        """
        return self._console_log_level

    def _create_file_handler(self) -> logging.FileHandler:
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setLevel(self.file_log_level)
        file_handler.setFormatter(logging.Formatter(
            fmt=self._file_fmt, datefmt=self._date_fmt))
        return file_handler

    def _create_console_handler(self, critical: bool) -> logging.StreamHandler:
        if critical:
            console_handler = logging.StreamHandler(sys.stderr)
            console_handler.addFilter(_CriticalFilter())
        else:
            console_handler = logging.StreamHandler(sys.stdout)
            console_handler.addFilter(_NonCriticalFilter())
        console_handler.setLevel(self.console_log_level)
        console_handler.setFormatter(logging.Formatter(fmt=self._console_fmt))
        return console_handler

    def get_logger(self, name: str) -> logging.Logger:
        r"""
        Get a logger by name
        """
        if name in self._loggers:
            return self._loggers[name]
        new_logger = logging.getLogger(name)
        new_logger.setLevel(logging.DEBUG)  # lowest level
        new_logger.addHandler(self._create_console_handler(True))
        new_logger.addHandler(self._create_console_handler(False))
        if self.log_file:
            new_logger.addHandler(self._create_file_handler())
        self._loggers[name] = new_logger
        return new_logger

    @log_file.setter
    def log_file(self, file_name: os.PathLike) -> None:
        self._log_file = file_name
        for logger in self._loggers.values():
            for idx, handler in enumerate(logger.handlers):
                if isinstance(handler, logging.FileHandler):
                    logger.handlers[idx].close()
                    if self.log_file:
                        logger.handlers[idx] = self._create_file_handler()
                    else:
                        del logger.handlers[idx]
                    break
            else:
                if file_name:
                    logger.addHandler(self._create_file_handler())

    @file_log_level.setter
    def file_log_level(self, log_level: int) -> None:
        self._file_log_level = log_level
        for logger in self._loggers.values():
            for handler in logger.handlers:
                if isinstance(handler, logging.FileHandler):
                    handler.setLevel(self.file_log_level)
                    break

    @console_log_level.setter
    def console_log_level(self, log_level: int) -> None:
        self._console_log_level = log_level
        for logger in self._loggers.values():
            for handler in logger.handlers:
                if type(handler) is logging.StreamHandler:  # pylint: disable=unidiomatic-typecheck
                    handler.setLevel(self.console_log_level)


log = LogManager()


def logged(obj: T) -> T:
    r"""
    Add logger as an attribute
    """
    obj.logger = log.get_logger(obj.__name__)
    return obj


#---------------------------- Interruption handling ----------------------------

@logged
class DelayedKeyboardInterrupt:  # pragma: no cover

    r"""
    Shield a code block from keyboard interruptions, delaying handling
    till the block is finished (adapted from
    `https://stackoverflow.com/a/21919644
    <https://stackoverflow.com/a/21919644>`__).
    """

    def __init__(self):
        self.signal_received = None
        self.old_handler = None

    def __enter__(self):
        self.signal_received = False
        self.old_handler = signal.signal(signal.SIGINT, self._handler)

    def _handler(self, sig, frame):
        self.signal_received = (sig, frame)
        self.logger.debug("SIGINT received, delaying KeyboardInterrupt...")

    def __exit__(self, exc_type, exc_val, exc_tb):
        signal.signal(signal.SIGINT, self.old_handler)
        if self.signal_received:
            self.old_handler(*self.signal_received)


#--------------------------- Constrained data frame ----------------------------

@logged
class ConstrainedDataFrame(pd.DataFrame):

    r"""
    Data frame with certain format constraints

    Note
    ----
    Format constraints are checked and maintained automatically.
    """

    def __init__(self, *args, **kwargs) -> None:
        df = pd.DataFrame(*args, **kwargs)
        df = self.rectify(df)
        self.verify(df)
        super().__init__(df)

    def __setitem__(self, key, value) -> None:
        super().__setitem__(key, value)
        self.verify(self)

    @property
    def _constructor(self) -> type:
        return type(self)

    @classmethod
    def rectify(cls, df: pd.DataFrame) -> pd.DataFrame:
        r"""
        Rectify data frame for format integrity

        Parameters
        ----------
        df
            Data frame to be rectified

        Returns
        -------
        rectified_df
            Rectified data frame
        """
        return df

    @classmethod
    def verify(cls, df: pd.DataFrame) -> None:
        r"""
        Verify data frame for format integrity

        Parameters
        ----------
        df
            Data frame to be verified
        """

    @property
    def df(self) -> pd.DataFrame:
        r"""
        Convert to regular data frame
        """
        return pd.DataFrame(self)

    def __repr__(self) -> str:
        r"""
        Note
        ----
        We need to explicitly call :func:`repr` on the regular data frame
        to bypass integrity verification, because when the terminal is
        too narrow, :mod:`pandas` would split the data frame internally,
        causing format verification to fail.
        """
        return repr(self.df)


#--------------------------- Other utility functions ---------------------------

def get_chained_attr(x: Any, attr: str) -> Any:
    r"""
    Get attribute from an object, with support for chained attribute names.

    Parameters
    ----------
    x
        Object to get attribute from
    attr
        Attribute name

    Returns
    -------
    attr_value
        Attribute value
    """
    for k in attr.split("."):
        if not hasattr(x, k):
            raise AttributeError(f"{attr} not found!")
        x = getattr(x, k)
    return x


def get_rs(x: RandomState = None) -> np.random.RandomState:
    r"""
    Get random state object

    Parameters
    ----------
    x
        Object that can be converted to a random state object

    Returns
    -------
    rs
        Random state object
    """
    if isinstance(x, int):
        return np.random.RandomState(x)
    if isinstance(x, np.random.RandomState):
        return x
    return np.random


@logged
def run_command(
        command: str,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        log_command: bool = True, print_output: bool = True,
        err_message: Optional[Mapping[int, str]] = None, **kwargs
) -> Optional[List[str]]:
    r"""
    Run an external command and get realtime output

    Parameters
    ----------
    command
        A string containing the command to be executed
    stdout
        Where to redirect stdout
    stderr
        Where to redirect stderr
    echo_command
        Whether to log the command being printed (log level is INFO)
    print_output
        Whether to print stdout of the command.
        If ``stdout`` is PIPE and ``print_output`` is set to False,
        the output will be returned as a list of output lines.
    err_message
        Look up dict of error message (indexed by error code)
    **kwargs
        Other keyword arguments to be passed to :class:`subprocess.Popen`

    Returns
    -------
    output_lines
        A list of output lines (only returned if ``stdout`` is PIPE
        and ``print_output`` is False)
    """
    if log_command:
        run_command.logger.info("Executing external command: %s", command)
    executable = command.split(" ")[0]
    with subprocess.Popen(command, stdout=stdout, stderr=stderr,
                          shell=True, **kwargs) as p:
        if stdout == subprocess.PIPE:
            prompt = f"{executable} ({p.pid}): "
            output_lines = []

            def _handle(line):
                line = line.strip().decode()
                if print_output:
                    print(prompt + line)
                else:
                    output_lines.append(line)

            while True:
                _handle(p.stdout.readline())
                ret = p.poll()
                if ret is not None:
                    # Handle output between last readlines and successful poll
                    for line in p.stdout.readlines():
                        _handle(line)
                    break
        else:
            output_lines = None
            ret = p.wait()
    if ret != 0:
        err_message = err_message or {}
        if ret in err_message:
            err_message = " " + err_message[ret]
        elif "__default__" in err_message:
            err_message = " " + err_message["__default__"]
        else:
            err_message = ""
        raise RuntimeError(
            f"{executable} exited with error code: {ret}.{err_message}")
    if stdout == subprocess.PIPE and not print_output:
        return output_lines