# Data

The datasets are public but too large for GitHub (~185k and ~370k cells), so they are not included in this repository.
Download them and place them as follows before running the notebooks:

```
data/
├── embryo/
│   ├── GSE157329_raw_counts.mtx.gz
│   ├── GSE157329_cell_annotate.txt.gz
│   └── GSE157329_gene_annotate.txt.gz
└── colorectal/
    ├── matrix.mtx.gz
    ├── genes.tsv
    └── barcodes.tsv
```

## 1. Human early embryogenesis (`notebooks/01_embryogenesis.ipynb`)

- Source: NCBI GEO, accession [GSE157329](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157329)
- Files: download the three supplementary files listed above (raw counts, cell annotation, gene annotation) and put them in `data/embryo/`.

## 2. Colorectal cancer immune hubs (`notebooks/02_colorectal_cancer.ipynb`)

- Source: Broad Institute Single Cell Portal, study [SCP1162](https://singlecell.broadinstitute.org/single_cell/study/SCP1162) (Pelka *et al.*, *Cell*, 2021)
- Files: download the raw count matrix (Matrix Market `.mtx.gz`) with its gene and barcode files, then rename them `matrix.mtx.gz`, `genes.tsv` and `barcodes.tsv` in `data/colorectal/`.
