# scVI virtual environment on Bayes
import numpy as np
import anndata as ad
import pandas as pd
import json
from scipy.io import mmwrite
import pandas as pd

adata = ad.read("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/yuhong/data/Gastrulation/erythroid_lineage.h5ad")
adata

#########

adata.obs.to_csv("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup14/pijunsala-erthyroid_adata_obs_backup.csv")
adata.obs = pd.DataFrame(index=adata.obs.index)  # Keep only the cell barcodes

#########

adata.var.to_csv("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup14/pijunsala-erthyroid_adata_var_backup.csv")
adata.var = pd.DataFrame(index=adata.var.index)  # Keep only gene name

#########

# Function to convert NumPy and other non-serializable objects into JSON-compatible formats
def convert_to_serializable(obj):
    if isinstance(obj, np.ndarray):  # Convert NumPy arrays to lists
        return obj.tolist()
    elif isinstance(obj, np.int64) or isinstance(obj, np.int32):  # Convert NumPy integers to Python integers
        return int(obj)
    elif isinstance(obj, np.float64) or isinstance(obj, np.float32):  # Convert NumPy floats to Python floats
        return float(obj)
    elif isinstance(obj, dict):  # Recursively handle dictionaries
        return {k: convert_to_serializable(v) for k, v in obj.items()}
    elif isinstance(obj, list):  # Recursively handle lists
        return [convert_to_serializable(v) for v in obj]
    elif isinstance(obj, tuple):  # Convert tuples to lists
        return list(obj)
    else:
        return obj  # Keep other types unchanged

# Convert AnnData uns metadata to a JSON-serializable format
uns_metadata_serializable = convert_to_serializable(adata.uns)

# Save to JSON file
output_path = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup14/pijunsala-erthyroid_adata_metadata.json"
with open(output_path, "w") as f:
    json.dump(uns_metadata_serializable, f, indent=4)

print(f"Metadata successfully exported to {output_path}")

adata.uns.clear()
adata

#########

# Ensure adata.raw exists
if adata.raw is not None:
    # Convert sparse matrix to Matrix Market (.mtx)
    mmwrite("/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup14/pijunsala-erthyroid_adata_raw_counts.mtx", adata.raw.X)
    print("Raw counts successfully exported as .mtx")
else:
    print("adata.raw is empty or does not exist.")

adata.raw = None  # Completely removes the raw attribute

########

# Directory to save the layer matrices
output_dir = "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup14/pijunsala-erthyroid_adata"

# Loop over all layers and export each as an .mtx file
if adata.layers is not None:
    for layer_name, layer_matrix in adata.layers.items():
        output_path = f"{output_dir}_layer_{layer_name}.mtx"
        mmwrite(output_path, layer_matrix)
        print(f"Layer '{layer_name}' successfully exported as .mtx")
else:
    print("adata.layers is empty or does not exist.")

adata.layers = None  # Completely removes the raw attribute

########

adata.write_h5ad(
    "/home/users/kzlin/kzlinlab/projects/veloUncertainty/out/kevin/Writeup14/pijunsala-erthyroid_cleaned.h5ad",
    compression='gzip', 
    compression_opts=9
)