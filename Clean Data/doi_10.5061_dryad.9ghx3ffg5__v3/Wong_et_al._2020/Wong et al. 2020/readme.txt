All data and R-scripts necessary to reproduce the results in Wong et al. (2020) are included here, except the phylogenetic tree of ants
The tree can be downloaded from the original paper (Economo et al. 2018) https://www.nature.com/articles/s41467-018-04218-4

A breif description of each file is provided here

.csv files
All_abundance.csv: Community data at different site. LMC = Lok Ma Chau, MP = Mai Po. Each cell (Column C - AE) represents the percentage of pitfall traps the species were captured.
                   However for fire ants (Solenopsis invicta) the number of workers was given instead. 
                   We are not aware of any null model analyses capable of dealing with percentage data, thus the code would convert the number to absence & presence in subsequent analyses.

env_data.csv: Environmental conditions at different sites. temp = temperature g.cover = bare ground cover

mean_trait.csv: mean trait values for all detected species. Size was log-transformed. mand = relative mandible length, pron = relative pronotum width, scap = relative scape length
                 
pdf_trait.csv: overlap of probability density function of each trait (i.e. absolute dissimilarity) of each species with S.invicta. 1 = no overlap; 0 = total overlap. The column pdf.all indicates overlap based on all traits.

.R files
Wong et al. 2020 main.R: all main functions required to produce the results

genus_tree2.R: randomly generate 100 species-level phylogenetic trees based on the genus-level tree from Economo et al. 2018

null_species_score_revised.R: obtain z-score by comparing null and observed networks

OR_ratio2.R: quantify strength of each co-occurrence relationships based on odd ratio (the source code is from the package spaa, but we modified it to conduct Haldane's correction.

species_summary2_revised.R: a function necessary for null_species_score_revised.R. It provides summary statistics for each species (e.g. total OR for all positive and negative links etc)

Please contact Toby Pak Nok Tsang (paknok@connect.hku.hk) if you have any questions! Thanks!!!