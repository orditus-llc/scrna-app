"""Backend for the Orditus guided Scanpy workflow.

Each module covers one analysis stage and contains computation only: no plotting,
printing or UI code. All parameter values live in `params.py`.

Every result is calculated from the input AnnData and the parameters passed in.
Expected PBMC3k results live only in `validation/` and `tests/`.
"""
