# Discovery sources (not trusted by the proof)

The C++ enumerator accepts `B I output_prefix`. The Python scripts archive the numerical discovery stage and require NumPy, SciPy, mpmath and (for metric14.py) numba. They are staging-directory tools, not part of `verify_all.py`. To rerun a stage, copy the relevant `core/` data/modules and `enumeration/` data into a separate working directory with these scripts, and decompress the enumeration files there. Do not run discovery in the proof directory: it may overwrite generated data.

The complete proof is replayed without any of these packages or discovery scripts. The finite data in `core/` and `enumeration/` are sufficient.
