# Generating Count Matrices with Snakemake

The CCA uses Snakemake to generate count matrices from FASTQ files. To process a new dataset, ensure you have the following inputs ready:

## Required Inputs

1. **FASTQ files**
   - Download FASTQ files to `data/dataset_id/fastqs`.

2. **Seqspec File**
   - Place the seqspec file in `data/dataset_id/fastqs`.
   - Include the specific "onlist" file for barcode error correction if necessary (refer to the [seqspec tutorial](https://github.com/pachterlab/seqspec) for details on generating a seqspec file).

3. **Observations File (`observations.txt`)**
   - The file should be tab-separated and include the tissue name and dataset ID.
   - The tissue must have a corresponding `markers.txt` file in the `markers/tissue` folder.

   Format: `dataset_id  tissue_name  technology_name  status`

4. **Snakemake Configuration**
- Add the dataset ID to the `RUNS` variable in the Snakemake file. e.g. `RUNS = ['dataset_id']`

## Required Programs

Ensure the following programs are installed (click on each program for installation instructions):
- [Snakemake](https://github.com/snakemake/snakemake)
- [gget](https://github.com/pachterlab/gget)
- [kb_python](https://github.com/pachterlab/kb_python)
- [jq](https://jqlang.github.io/jq/download/)
- [mx](https://github.com/cellatlas/mx)
- [ec](https://github.com/cellatlas/ec)
- [seqspec](https://github.com/pachterlab/seqspec)

## Generating Matrices

Run the following command to generate the count matrices:

```bash
snakemake -c 64 # -c sets the number of threads
```
