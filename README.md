# Correlation-misspecification-in-longitudinal-gene-expression-studies
 This GitHub page serves as a repository for the results provided in the manuscript titled "Correlation misspecification in longitudinal gene expression studies" and its supplementary text.
The repository is structured as a sequence of R scripts each labeled by the corresponding figure or table found in the manuscript.  We verified that the codes run successfully on R version 4.2.0.  The figures and tables that each R script produces are provided in the "Figs and Tables" folder for easy access.  Additional inquiries into the codes and outputs can be made by contacting the corresponding author at turnerja2@sfasu.edu.   

Note 1: It should be noted that each of the figures that provide simulation results can be time-consuming. In particular, the FDR simulations are provided in Figure 3 and Figure 6.  Reviewers might find it helpful to do a fewer number of simulations first to verify the code is working as intended.

Note 2: We do not provide the raw data files for the Flu challenge data set used as an illustration of our method in the manuscript.  Because the files are quite large, we provided the necessary coding steps to extract the data from the Gene Expression Omnibus using Bioconductor.  Both the expression files and the experimental design file are then imported into data frames, cleaned, and processed for easily modeling of each transcript (gene) one at a time.  These steps can be found at the beginning of the Rscript for Figures 7, 8, and Table 1.
