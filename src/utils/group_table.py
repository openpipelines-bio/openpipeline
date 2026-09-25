"""Helpers for group-level (e.g. donor / sample) tables.

BEYOND-style analyses operate on a group x label matrix (participants x
subpopulation proportions), not on cells. Such a matrix does not belong in a
MuData: the cell-level data is not used at all. These helpers read and write it
as a plain CSV whose first column holds the group identifier.
"""

import pandas as pd


def read_group_table(path, id_column=None, label="Input table"):
    """Read a group x feature CSV.

    The identifier column is `id_column`, or the first column when `id_column`
    is None. Every remaining column must be numeric.

    Returns a DataFrame indexed by the group identifier, with the name of the
    identifier column kept as `index.name`.
    """
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"{label} '{path}' has no rows.")
    if df.shape[1] < 2:
        raise ValueError(
            f"{label} '{path}' has {df.shape[1]} column(s); expected an "
            "identifier column plus at least one value column."
        )

    if id_column is None:
        id_column = df.columns[0]
    elif id_column not in df.columns:
        raise ValueError(
            f"Identifier column '{id_column}' not found in {label} '{path}'. "
            f"Available: {list(df.columns)}"
        )

    ids = df[id_column].astype(str)
    if ids.duplicated().any():
        duplicated = sorted(ids[ids.duplicated()].unique())
        raise ValueError(
            f"{label} '{path}' has duplicated identifiers in column "
            f"'{id_column}': {duplicated}"
        )

    values = df.drop(columns=[id_column])
    non_numeric = [
        col for col in values.columns if not pd.api.types.is_numeric_dtype(values[col])
    ]
    if non_numeric:
        raise ValueError(
            f"{label} '{path}' has non-numeric value column(s): {non_numeric}. "
            f"Only '{id_column}' may be non-numeric."
        )

    values.index = pd.Index(ids.values, name=id_column)
    return values


def write_group_table(df, path, id_column=None):
    """Write a DataFrame indexed by group identifier to CSV.

    The index becomes the first column, named `id_column` or `df.index.name`,
    falling back to `group`.
    """
    out = df.copy()
    name = id_column or out.index.name or "group"
    out.index = pd.Index(out.index, name=name)
    out.to_csv(path, index=True)
    return name
