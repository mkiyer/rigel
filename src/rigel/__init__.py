from importlib.metadata import version as _version, PackageNotFoundError

try:
    __version__ = _version("hulkrna")
except PackageNotFoundError:
    __version__ = "0.0.0+unknown"
