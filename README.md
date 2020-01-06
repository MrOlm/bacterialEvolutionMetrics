# bacterialEvolutionMetrics

Scripts in support of the publication "Consistent Metagenome-Derived Metrics Verify and Delineate Bacterial Species Boundaries"

Open source pre-print publication is available at
[bioRxiv](https://www.biorxiv.org/content/10.1101/647511v1)

Publication is in press at mSystems

An open-source program enabling marker gene analysis from metagenomic data is available at https://github.com/alexcritschristoph/RPxSuite

Nucleotide sequences of genome sets used in this study are available at https://doi.org/10.6084/m9.figshare.c.4508162.v1

## Reproducing figures in publication

In the JupyterNotebooks folder is a list of Jupiter Notebooks. These can be viewed using [Jupiter](https://jupyter.org/), and have the details on how much of the analysis was performed.

Notebook 1 documents how to load the results of FastANI in a standardized manner.

Notebook 2 documents how to use the output from Notebook 1 to generate Figure 1

Notebook 3 documents how to load the results of dnds_from_drep.py, which is needed to make Figure 2. dnds_from_drep.py is a python script that is to be run on the output from from dRep when run using the `--S_algorithm goANI` option.

Notebook 4 documents how to run popCOGent in a hack-y way. Methods are largely copy pasted from the original code [here](https://github.com/philarevalo/PopCOGenT/tree/master/src/PopCOGenT), based on the following [publication](https://linkinghub.elsevier.com/retrieve/pii/S0092867419307366). This notebook allows you to run only the part of the program that is needed to generate Figure 2

Notebook 5 takes the results of Notebooks 3 and 4 to generate Figure 2.

Notebooks 6 and 7 document how species delineation thresholds and Figure 3 were generated.
