# SMAC, the Smart, Automatic Classification system

## Introduction
The **SM**art, **A**utomatic **C**lassification system (SMAC) aims at integrating and analysing literature, biomedical and experimental data starting from a user-defined, text-free search query

SMAC exploits the [Entrez Programming Utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/) to mine biomedical literature and identify the most relevant articles. Records are then ranked according to the [“Best Match” relevance algorithm](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.2005343).

Medical Subject Headings (MeSH) are ranked according to their (i) hierarchical level (specificity), (ii) abundance in topics-related articles and (iii) abundance in all PubMed citations. Firstly, MeSH terms are sorted separately according to each criteria, then the [Robust Ranking Aggregation method](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3278763/) is applied.

## Starting the analysis

SMAC has been developed using Python and R, and is distributed to public in a form of [Docker package](https://hub.docker.com/r/hfx320/smac). Thanks to its modularity, users can easily edit the source code of SMAC by implementing R-based, custom analyses to be included in the main pipeline.

### Pulling the docker image from Docker Hub
**Please note. A working network connection is required**
To download the latest SMAC docker, use the command `docker pull hfx320/smac`. After few minutes (depending on your network speed), SMAC will be automatically installed and configured on your machine, ready to be used.

### Running SMAC
Once SMAC docker is correctly downloaded and installed on your system, the following command can be used to run your first analysis:

`docker run -v <folder where to save the results of your analysis>:/smac/results hfx320/smac:latest -t "breast AND CANCER" -l 10000 -s 1 -e <your email address>`

By running this command (once defined a valid folder and email address), SMAC will explore the literature and will retrieve the first impacting 10,000 papers about BREAST CANCER. Since we skipped the retrieval of molecular data (`-s 1`), the analysis will not continue further.

A detailed description of all the available commands for SMAC is available while typing:
`docker run -v <folder where to save the results of your analysis>:/smac/results hfx320/smac:latest -h`

`-t TERMS [TERMS ...], --terms TERMS [TERMS ...]
                        terms for performing the analysis
  -p PMIDS, --pmids PMIDS
                        file with list of pmids to download and analyse
  -l LIMIT, --limit LIMIT
                        max number of publications to retrieve
  -a ANALYSIS, --analysis ANALYSIS
                        analysis to perform (separated by comma) [pca,receptor
                        _status,molecular_classification,tumour_purity,gene_ex
                        pression]
  -s {0,1}, --skip {0,1}
                        do you want to skip GSE download/analysis [1=yes, skip
                        the analysis, 0=no]
  -e EMAIL, --email EMAIL
                        email for retrieving the papers
`
