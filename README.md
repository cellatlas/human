# Human Commons Cell atlas
This GitHub reposotiry for the Human Commons Cell Atlas contains uniformly preprocessed and filtered count matrices and cell type assignments that span 27 Human organs as well as marker gene lists for 31 Human organs. The count matrices are in the `data` folder and the marker gene lists in the `markers` folder. 

The repository has the following structure (note: the repository is large):

```bash
├── data
│   ├── adipose
│   │   ├── GSM3711757
│   │   │   ├── metadata.json
│   │   │   ├── preprocess.ipynb
│   │   │   ├── filter
│   │   │   │   ├── barcodes.txt.gz
│   │   │   │   ├── genes.txt.gz
│   │   │   │   └── matrix.mtx.gz
│   │   │   ├── metrics
│   │   │   │   ├── inspect.json
│   │   │   │   ├── kb_info.json
│   │   │   │   ├── kneeplot.png
│   │   │   │   ├── mx_metrics.json
│   │   │   │   └── run_info.json
│   │   │   └── mx_out
│   │   │       └── assignments_rank_mx.tsv
│   │   ├── GSM3711758
│   │   ├── GSM3711759
|   |   ...
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
└── human_datasets.csv
```

Summary statistics for the Human Commons Cell Atlas can be found here: https://cellatlas.github.io/human
