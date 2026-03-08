import os as _os

# Prevent numpy/OpenMP from spawning idle thread pools that compete with
# rigel's own C++ parallelism for cache and scheduling resources.
_os.environ.setdefault("OMP_NUM_THREADS", "1")

from importlib.metadata import version as _version, PackageNotFoundError

try:
    __version__ = _version("rigel-rnaseq")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"
