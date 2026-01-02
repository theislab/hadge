#!/usr/bin/env python3

# versions
import platform
import yaml

import pandas as pd
import scanpy as sc
import mudata as md
import anndata as ad

from anndata import AnnData
from pathlib import Path
from mudata import MuData


class Arguments:
    """Parses the arguments, including the ones coming from $task.ext.args.
    Adopted from mygene module (Suzanne Jin)."""

    def __init__(self) -> None:
        self.parse_input_args()
        self.creat_output_dirs()

    def parse_input_args(self) -> None:
        self.prefix = "$task.ext.prefix" if "$task.ext.prefix" != "null" else "$meta.id"
        self.negative_str = "negative"

        self.hto_matrix = "${hto_matrix}"
        self.rna_matrix = "${rna_matrix}"
        self.hashing_summary_assignment = "${hashing_summary_assignment}"
        self.hashing_summary_classification = "${hashing_summary_classification}"
        self.genetic_summary_assignment = "${genetic_summary_assignment}"
        self.genetic_summary_classification = "${genetic_summary_classification}"

        path_vars = {
            "rna_matrix",
            "hto_matrix",
            "hashing_summary_assignment",
            "hashing_summary_classification",
            "genetic_summary_assignment",
            "genetic_summary_classification",
        }

        def _tranlate_to_python(input_str, value_str):
            if value_str.strip() == "":
                return None
            else:
                path = Path(value_str)
                if not path.exists():
                    raise FileNotFoundError(f"Path does not exist: {path}")
                return path

        for var in path_vars:
            raw_value = getattr(self, var)
            processed_value = _tranlate_to_python(var, raw_value)
            setattr(self, var, processed_value)

    def creat_output_dirs(self) -> None:
        directories = {
            "h5mu": "_genetic_and_hashing.h5mu",
            "genetic": "_genetic.h5ad",
            "hashing": "_hashing.h5ad",
        }

        for output, directory in directories.items():
            setattr(self, output, self.prefix + directory)

    def print_args(self) -> None:
        for attr in vars(self):
            print(f"{attr}: {getattr(self, attr)}")


def saveAnnData(args: Arguments, isRNA: bool, count_data: AnnData):
    if isRNA:
        summary_files = {
            "genetic_summary_assignment",
            "genetic_summary_classification",
        }
    else:
        summary_files = {
            "hashing_summary_assignment",
            "hashing_summary_classification",
        }

    for file in summary_files:
        path = getattr(args, file)
        if path is None:
            continue
        summary_table = pd.read_csv(path)
        summary_table.set_index("Barcode", inplace=True)
        summary_table = summary_table.add_suffix(f"_{file.split('_')[-1]}")
        count_data.obs = count_data.obs.join(summary_table, how="left").fillna(
            args.negative_str
        )

    if isRNA:
        count_data.write(args.genetic)
    else:
        count_data.write(args.hashing)

    return count_data


if __name__ == "__main__":
    args = Arguments()

    if args.rna_matrix is not None and args.hto_matrix is not None:
        rna_data = sc.read_10x_mtx(args.rna_matrix)
        hto_data = sc.read_10x_mtx(args.hto_matrix, gex_only=False)

        rna_data = saveAnnData(args, True, rna_data)
        hto_data = saveAnnData(args, False, hto_data)

        mudata = MuData({"rna": rna_data, "hto": hto_data})
        mudata.update()
        mudata.write(args.h5mu)

    elif args.rna_matrix is not None:
        rna_data = sc.read_10x_mtx(args.rna_matrix)
        saveAnnData(args, True, rna_data)

    elif args.hto_matrix is not None:
        hto_data = sc.read_10x_mtx(args.hto_matrix, gex_only=False)
        saveAnnData(args, False, hto_data)

    versions = {
        "${task.process}": {
            "python": platform.python_version(),
            "pandas": pd.__version__,
            "scanpy": sc.__version__,
            "mudata": md.__version__,
            "anndata": ad.__version__,
        }
    }

    with open("versions.yml", "w") as f:
        yaml.dump(versions, f)
