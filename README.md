# COX_Proteome-stability

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4280621.svg)](https://doi.org/10.5281/zenodo.4280621)


This repository contains the analysis code associated with the Proteome Folding Stability project, led by Dr. Dezerae Cox. This manuscript has been submitted for publication under the title "Surveyance of proteome foldedness unmasks hidden information on protein binding and conformation".

Final version will be provided upon publication.

## Prerequisites

This analysis assumes a standard installation of Python 3 (=> 3.6). For specific package requirements, see the environment.yml file, or  create a new conda environment containing all packages by running ```conda create -f environment.yml```. In addition to the analysis contained here, some simple statistical tests were performed using [GraphPad Prism v 8.0](https://www.graphpad.com/scientific-software/prism/). 

## Raw data and databases

In the case of the published dataset comparison, raw data in the form of supplementary datasets from each publication can be accessed in part using the appropriate ```raw_data.py``` script. Unfortunately, due to journal subscription requirements, in some instances journal access is required. In these cases, datasets will need to be manually downloaded.

Initial processing of the novel mass spectrometry spectra files was completed using either Proteome Discoverer or MaxQuant. The preprocessed identification and quantitation data have been deposited alongside the ```.RAW``` files via the PRIDE [1] partner repository to the ProteomeXchange Consortium under the dataset identifiers PXD022587 and PXD022640. For convenience, the preprocessed identification and quantitation data (hereon termed raw data) have also been uploaded alongside the other raw data as an open-access [Zenodo dataset](https://doi.org/10.5281/zenodo.4280620). These datasets can be collected automatically using the ```raw_data.py``` script in each of the respective analysis folders.

Lastly, various public databases were queried as indicated in the accompanying manuscript, for which access protocols are also provided in the respective analysis workflow where appropriate.

## Workflow

To reproduce analyses presented in the manuscript, it is recommended to use the following order: published dataset comparison, lysate denaturation, recombinant client assay. This will ensure any source data is generated as needed (e.g. lysate denaturation requires summary of the published datasets to measure the correlation). In addition, where processing order is important for individual analyses, scripts have been numbered and should be run in order before unnumbered counterparts.

## A quick note on cluster numbers

By default, python is zero-indexed and therefore the clusters are labelled 0-3 initially (rather than 1-4, as in the manuscript). In addition, the assigned number during the initial clustering is stochastic (as it is unsupervised, the cluster centres can be initialised in a randomized order). However, the peptides associated to each cluster should be conserved. Therefore, cluster numbers were manually mapped to the order presented in the manuscript within the analysis workflow (cluster map dictionary in ```cluster_summary.py``` script).


## References

[1]. Perez-Riverol Y, Csordas A, Bai J, Bernal-Llinares M, Hewapathirana S, Kundu DJ, Inuganti A, Griss J, Mayer G, Eisenacher M, Pérez E, Uszkoreit J, Pfeuffer J, Sachsenberg T, Yilmaz S, Tiwary S, Cox J, Audain E, Walzer M, Jarnuczak AF, Ternent T, Brazma A, Vizcaíno JA (2019). The PRIDE database and related tools and resources in 2019: improving support for quantification data. Nucleic Acids Res 47(D1):D442-D450 (PubMed ID: 30395289)