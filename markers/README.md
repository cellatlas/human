# Markers 📍

Google Colab notebooks download and process organ-specific cell-type markers from published literature.

## Adding markers

The CCA atlas is a "live" atlas that can be iteratively built and improved with updated cell type annotations. Markers are derived from tables of differentially expressed genes provided in published literature. To add new markers to the atlas from published literature, identify the following information from a paper:

- `species`: the latin name of the species, lowercase, separated by underscores (e.g. homo_sapiens)
- `organ`: the name of the organ, lowercase, separated by underscores (e.g. mammary_gland)
- `reference`: the specific reference used to generated the differential expression table (e.g. GRCh37-Enesmbl75)
- `paper_doi`: the doi of the paper (e.g. https://doi.org/10.1016/j.devcel.2020.05.010)
- `table_link`: the direct link to the differential expression table (e.g. https://ars.els-cdn.com/content/image/1-s2.0-S1534580720303993-mmc2.xlsx)

To produce cell-type markers from that paper (e.g. `markers.txt`) create a Google Colab notebook that performs the following tasks:

1. Downloads the table
2. Filter genes by certain criteria (this depends on the columns of the table- e.g. pvalue, logfc, rank)
3. Group genes by celltype
4. Write `markers.txt` to disk

Importantly, declare define the following object in the notebook so that it can be included in the header of the `markers.txt` file.

```python
header = [
    {
      "species": species,
      "organ": organ,
      "reference": reference,
      "paper_doi": paper_doi,
      "table_link": table_link,
    }
]
```

Note: if you are pulling from multiple sources, consider either creating a new notebook for the second source, or adding another entry to the `header` object:

```python
header = [
    {
      "species": species,
      "organ": organ,
      "reference": reference_1,
      "paper_doi": paper_doi,
      "table_link": table_link_1,
    },
    {
      "species": species,
      "organ": organ,
      "reference": reference_2,
      "paper_doi": paper_doi,
      "table_link": table_link_2,
    }
]
```

For an example of how to do generate a `markers.txt` file, see this Google Colab notebook: [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cellatlas/human/blob/master/markers/testis/markers.ipynb)
