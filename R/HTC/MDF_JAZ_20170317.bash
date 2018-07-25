#!/bin/csh
#$ -M jzwart@nd.edu
#$ -m abe
#$ -pe mpi-24 48
#$ -q long
#$ -N MDF_JAZ_20170317
#$ -r y

module load bio/R

setenv R_LIBS ~/NHLDLakeCarbonModel/R_packages

Rscript fakeData_2DOCpools_MM_Temp_20170317.R
