#!/usr/bin/env python3
import scanpy as sc

input_mtx_dir = "${input_mtx_dir}"
prefix = "${prefix}"

# Read 10X matrix
adata = sc.read_10x_mtx(input_mtx_dir,gex_only=False)

# Save as h5ad
output_path = f"{prefix}_hto_data.h5ad"
adata.write_h5ad(output_path)

# Version info
with open("versions.yml", "w") as f:
    f.write("this is a test")
