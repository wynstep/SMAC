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

### Available types of analysis
The updated version of SMAC provides a wide range of ways for exploring the scientific literature:

* __MeSH term__ : This is the classical way for conducting an analysis, the input will be a sentence composed by terms, bound by logic gates (AND, OR, NOT). For example, if I want to retrieve papers and data associated to "cancer and microRNAs", I would just have to type: __cancer AND microRNAs__.<br>
__PLEASE NOTE__ : SMAC can automatically extract the synonyms for a selected term (case insensitive). As consequence, in this specific example we don't have to declare synonyms of microRNAs (i.e. _miRNAs_) or cancer (i.e. _tumour_).

* __PMID__: If the user already has a list of IDs from PubMed and wants to download molecular data as well as perform bioinformatics analyses on these, SMAC has made this type of analysis available. In this context, the input of our analysis will be a text file (any extension accepted) containing a list of PMIDs to parse.

* __GeneName__: In case the user is interested in extracting the most relevant informations related to a single/list of genes, it's now possible to perform this analysis using SMAC.
As for the PMID type of analysis, the input is a text file (any extension accepted) that contains a list of gene names (HGNC format).

### Available options
SMAC provides a wide range of options, according to the selected type of analysis.

##### MeSH Terms
`docker run -v <folder where to save the results of your analysis>:/SMAC2/DATA hfx320/smac:latest -t term -p "-h"`
| Short Flag | Long Flag | Values | Description | Required |
| ------------ | ------------- | ------ | ----------- | ----- |
| -t           | --queryTerm       | Any sentence separated by AND/OR logic gates | Term for starting a PubMed search query | TRUE |
| -l | --PMIDsLimit | Int number | Max number of PMIDs to retrieve | FALSE <default = all> |
| -a | --analysis | String | analysis to perform (separated by comma) [pca,receptor_status,molecular_classification,tumour_purity,gene_expression,functional_enrichment] | FALSE <default = pca,gene_expression> |
| -s | --skipDownload | String | Do you want to skip download of molecular data? <yes,no> | FALSE <default = no> |
| -g | --nGenes | Int number | Number of top differentially-expressed genes to plot in the gene expression analysis | FALSE <default = 100> |
| -n | --nCores | Int number | Number of cores/threads for the analysis | FALSE <default = 1> |
| -e | --email | String | A valid email address for not being banned by NCBI | TRUE |

##### PubMed IDs
`docker run -v <folder where to save the results of your analysis>:/SMAC2/DATA hfx320/smac:latest -t pmid -p "-h"`
| Short Flag | Long Flag | Values | Description | Required |
| ------------ | ------------- | ------ | ----------- | ----- |
| -p | --PMIDlist | File | List of PMIDs to parse | TRUE |
| -a | --analysis | String | analysis to perform (separated by comma) [pca,receptor_status,molecular_classification,tumour_purity,gene_expression,functional_enrichment] | FALSE <default = pca,gene_expression> |
| -s | --skipDownload | String | Do you want to skip download of molecular data? <yes,no> | FALSE <default = no> |
| -g | --nGenes | Int number | Number of top differentially-expressed genes to plot in the gene expression analysis | FALSE <default = 100> |
| -n | --nCores | Int number | Number of cores/threads for the analysis | FALSE <default = 1> |
| -e | --email | String | A valid email address for not being banned by NCBI | TRUE |

##### Gene names
`docker run -v <folder where to save the results of your analysis>:/SMAC2/DATA hfx320/smac:latest -t gene -p "-h"`
| Short Flag | Long Flag | Values | Description | Required |
| ------------ | ------------- | ------ | ----------- | ----- |
| -t | --geneList | File | List of PMIDs to parse | TRUE |
| -l | --PMIDsLimit | Int number | Max number of PMIDs to retrieve | FALSE <default = all> |
| -a | --analysis | String | analysis to perform (separated by comma) [pca,receptor_status,molecular_classification,tumour_purity,gene_expression,functional_enrichment] | FALSE <default = pca,gene_expression> |
| -s | --skipDownload | String | Do you want to skip download of molecular data? <yes,no> | FALSE <default = no> |
| -g | --nGenes | Int number | Number of top differentially-expressed genes to plot in the gene expression analysis | FALSE <default = 100> |
| -n | --nCores | Int number | Number of cores/threads for the analysis | FALSE <default = 1> |
| -e | --email | String | A valid email address for not being banned by NCBI | TRUE |

### Running SMAC
Once SMAC docker is correctly downloaded and installed on your system, you're ready to go! <br>
This user guide will cover some examples to run fast analyses and test SMAC functionalities

##### Retrieve data associated to BREAST CANCER
`docker run -v <folder where to save the results of your analysis>:/SMAC2/DATA hfx320/smac:latest -t term -p "--queryTerm 'breast AND cancer' -l 10000 -s 1 -e <email address>"`

By running this command (once defined a valid folder and email address), SMAC will explore the literature and will retrieve the first impacting 10,000 papers about BREAST CANCER. Since we skipped the retrieval of molecular data (`-s 1`), the analysis will not continue further.

##### Retrieve data associated to TP53, KRAS and SMAD4 genes
`docker run -v <folder where to save the results of your analysis>:/SMAC2/DATA hfx320/smac:latest -t gene -p "--geneList <txt file containing genes> -l 1000 -s 0 -g 50 -a pca,gene_expression,functional_enrichment -n 4 -e <email address>"`

This command will download the top 1,000 articles and relative molecular data (if available) that are connected with TP53, SMAD4 and KRAS. The bioinformatics analyses conducted on the downloaded data will be __PCA, differential gene expression (top 50 genes up/down regulated), functional enrichment__. The process will be performed on 4 threads.
