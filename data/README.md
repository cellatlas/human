# Observations 👀

Data are organized by organ. A single-cell count matrix derived from a set of FASTQ files is labeled an observation. Multiple observations have been compiled for each organ.

Each observation contains a filtered gene count matrix, cell type labels and metadata as well as code to produce these files.

## Files

1. `barcodes.txt.gz`: list of barcodes (size n)
2. `genes.txt.gz`: list of gene names (size m)
3. `gene_ids.txt.gz`: list of gene ids (size m)
4. `matrix.mtx.gz`: a barcode by genes matrix (size n x m)
5. `labels.txt.gz`: barcode to cell type map (provided by paper)
6. `metadata.json`: metadata for the observation
7. `kb_info.json`: metadata for the whole pre-processing workflow
8. `inspect.json`: metadata for the BUS file
9. `run_info.json`: metadata for pseudoalignment
10. `assignments.txt.gz`: cell type assignments (size n)
11. `groups.txt`: cell types present in dataset

## Code
- `preprocess.ipynb`: notebook to produce matrices and metadata
  - In:
    - SRA / ENA Accession
  - Out:
    1. `barcodes.txt.gz`
    2. `genes.txt.gz`
    3. `gene_ids.txt.gz`
    4. `matrix.mtx.gz`
    5. `labels.txt`
    6. `metadata.json`
    7. `kb_info.json`
    8. `inspect.json`
    9. `run_info.json`
- `assign.ipynb`: notebook to produce cell assignements based on markers
  - In:
    1. `barcodes.txt.gz`
    2. `genes.txt.gz`
    3. `gene_ids.txt.gz`
    4. `matrix.mtx.gz`
    5. `markers.txt` (from `../markers/<ORGAN>/`)
  - Out:
    1. `assignments.txt.gz`
    2. `groups.txt`
    
