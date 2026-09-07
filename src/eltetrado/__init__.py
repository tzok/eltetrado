try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:
    from importlib_metadata import PackageNotFoundError, version

try:
    __version__ = version("eltetrado")
except PackageNotFoundError:
    __version__ = "0.0.0.dev0"
