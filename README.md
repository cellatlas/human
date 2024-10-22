# Human Commons Cell Atlas
This GitHub repository for the Human Commons Cell Atlas contains uniformly preprocessed and filtered count matrices and cell type assignments that span 29 Human organs as well as marker gene lists for 31 Human organs. The count matrices are in the `data` folder and the marker gene lists in the `markers` folder. 

All matrices where aligned to the [Ensembl v96 reference](https://github.com/pachterlab/kallisto-transcriptome-indices/releases/tag/ensembl-96) which corresponds to Gencode v30 and UCSC v30 and Genome Assembly version GRCh38.p12

![alt text](https://github.com/cellatlas/human/blob/main/docs/vetruvian_man.png?raw=true)

The repository has the following structure (note: the repository is large):

```bash
├── data
│   ├── GSM3711757
│   │   ├── metadata.json
│   │   ├── preprocess.ipynb
│   │   ├── kb_out_nac
│   │   └── mx_out
│   │       └── assignments.txt
│   ├── GSM3711758
│   ├── GSM3711759
|   ...
├── tissues
│   ├── adipose
│   │   ├── GSM3711757 # symlink to data folder
│   │   ├── ...
│   ├── bladder
│   ├── blood
│   ├── bone_marrow
│   ├── brain
|   ...
├── markers
│   ├── adipose
│   │   ├── markers.ipynb
│   │   └── markers.txt
│   ├── bladder
│   ├── blood
│   ├── bone_marrow
│   ├── brain
|   ...
├── docs
│   ├── index.html
│   └── metrics.json
├── observations.txt # (to be combined with datasets.tsv soon)
└── datasets.tsv # (to be combined with observations.txt soon)
```

Summary statistics for the Human Commons Cell Atlas can be found here: https://cellatlas.github.io/human
