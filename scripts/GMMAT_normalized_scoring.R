
#####################################
###Load in data and set variables####
#####################################

#set wd, load all necassary libraries
setwd("/home/bal7cg/R/data.objects/")
getwd()

library(data.table)
library(SeqArray)
library(GMMAT)
library(stringr)
library(tidyverse)
library(MAGEE)
#Bring in command line arguments to define job ID
args <- commandArgs(trailingOnly=TRUE)#load in expression data
#jobID <- as.numeric(args[1]) 
jobID <- 1

#load in GRM
GRM <- readRDS("GRM_All")
#load in expression data(mean reads > 2), normalized with edgeR trimmed means of M values
expression_matrix_filename <- readRDS("tmmcounts")


#call the gene name for each jobID
geneID <- expression_matrix_filename[jobID,1]
#create filenames for null models and score report
null_model_filename <- paste("/scratch/bal7cg/models/","model",geneID, sep ="")
score_report_output <- paste("/scratch/bal7cg/score_output/","score",geneID, sep ="")

gds_filename <- "/home/bal7cg/R/data.objects/tmp.gds"

showfile.gds(closeall = TRUE)



#####################################################
#create pheno file for a gene based on Job ID####
####################################################
#seperate out gene of interest

Gene_int_exp <- expression_matrix_filename[jobID,]

#create a data.table with expression data, pre/post treatment, and cage numbers listed for every sample
Genemelt= melt(Gene_int_exp,id.vars = c("GeneID"),variable.name= "samplename", value.name = "reads")


Genemelt[, c("pre_pos","cage","line") := tstrsplit(samplename,"[.]")]
#Turn pre/post names into 0/1 integer values
Treatment <- str_replace_all(Genemelt$pre_pos,"4CB","0")
Treatment <- str_replace_all(Treatment,"23HR","1")
Treatment <- as.integer(Treatment)

#change samplename periods to hyphens (to match gds)
newnames <- Genemelt[[2]]
newnames <- c(gsub(".","-",newnames, fixed = T))

#change cage numbers to integers (for glmmkin)
integer_cage <- as.integer(Genemelt$cage)
#reorganize data.table
pheno <- data.table(
  Genemelt,
  Treatment,
  newnames,
  integer_cage
)
pheno <- pheno[,c(1,8,9,7,3)]


####################
##Create Null Model#
####################

#fit a GLMM to the data
null_model <- glmmkin(reads~integer_cage+Treatment, data = pheno, kins = GRM, id = "newnames",family = poisson(link = "log"))


###########################
###Run a Score Report######
##########################
#run a Gene by environment score report
glmm.gei(null_model, interaction = c("Treatment", "integer_cage"), geno.file = gds_filename,outfile = score_report_output,ncores = 10)
#Run a GMMAT score report using the null_model and gds file


