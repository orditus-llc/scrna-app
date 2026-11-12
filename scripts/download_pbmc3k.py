from pathlib import Path
import scanpy as sc

# Find the main project folder
project_folder = Path(__file__).resolve().parent.parent

# Create the path to the data folder
data_folder = project_folder / "data"

# Download the raw PBMC3k dataset
adata = sc.datasets.pbmc3k()

# See what we downloaded
print(adata)

# Save it as an AnnData (.h5ad) file
output_file = data_folder / "pbmc3k.h5ad"
adata.write_h5ad(output_file)

print(f"Saved PBMC3k to: {output_file}")