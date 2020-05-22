#Iterated_scoring Gather
#####################################
##load in data and set variables####
###################################
setwd("/home/bal7cg/R/data.objects/")
library(data.table)
library(tidyverse)

library(foreach) ### AOB
library(doMC) ### AOB
library(heritability)
registerDoMC(10)



#bring in the normalized expression data set, use it as key to call score reports
expression_matrix_filename <- readRDS("tmmcounts")

####################################
##construct heritability data table##
####################################
#read in the already gathered SNP data
#split_out <- readRDS("split_out")
#read in GRM that has same sample names as expression data
GRM <- readRDS("GRM_All")

#filter out the gene expression data to only include genes from the first gather
#short_genes <- split_out %>% 
 # distinct(geneid)
#generows <- match(short_genes$geneid,expression_matrix_filename$GeneID)

#filtered_expression <- expression_matrix_filename[generows,]
filtered_genes <- expression_matrix_filename$GeneID
#shortgenes <- filtered_genes[1:3]
####################################
##Use foreach loop for each gene  ##
####################################
#create a foreach loop that runs through each gene
#uses parallel computing accross 10 cores

out <- foreach(f=filtered_genes)%dopar%{
  message(f)
  
  ################################ 
  ##create expression data table##
  ################################
  #This part of the loop creates a table with expression as phenotype, 
  #and includes factors such as cage, pre/pos treatment as covariates
  #same code as was used for initial Slurm Job
  
  #seperate out gene of interest
  
Gene_int_exp <- expression_matrix_filename %>% 
    filter(GeneID == "FBgn0030004")
  
  #create a data.table with expression data, pre/post treatment, and cage numbers listed for every sample
  Gene_int_exp <- as.data.table(Gene_int_exp)
  Genemelt= melt(Gene_int_exp,id.vars = c("GeneID"),variable.name= "samplename", value.name = "reads")
  Genemelt = Genemelt[, c("pre_pos","cage","line") := tstrsplit(samplename,"[.]")]
  
  #Turn pre/post names into 0/1 integer valuessd
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
  pheno <- pheno[,c(1,2,8,9,7,3)]
  
  ########################
  ##Measure heritability##
  ########################
  #covariates <- pheno[,c(4,5)]
  cage <- pheno[[4]]
  treatment <- pheno[[5]]
  covariates <- cbind(cage,treatment)
  heritable <- marker_h2(data.vector = pheno[[6]], geno.vector = pheno[[3]], covariates = pheno[[5]],K = GRM, max.iter = 10)
  
  x <- data.frame("geneid" = f, "va" = heritable$va, "ve" = heritable$ve, "h2" = heritable$h2)
}
#turn list into data.table
normalized.heritability <- rbindlist(out)

saveRDS(normalized.heritability, "normalized.heritability")
norm <- readRDS("normalized.heritability")
#compare pre-normalized heritability scores with post normalized
norm <- norm[,normal := "normalized"]
nonorm <- readRDS("heritability_out")
nonorm <- nonorm[,normal := "notnormalized"]
heridata <- rbind(nonorm,norm)
heridata <- transform(heridata, normal = as.factor(normal))
ggplot(heridata, aes(geneid,h2, color = normal)) + geom_point(alpha = .6, shape = ".", size = 3)
ggplot(heridata, aes(x = normal, y = h2)) + 
  geom_boxplot(fill = "darkmagenta", color = "#C4961A")
#some of the normalized scores are unusually high.
normal[h2 >.6]
heridata[geneid == "FBgn0030004"]
heridata[geneid == "FBgn0038371"]
