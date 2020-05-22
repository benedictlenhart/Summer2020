###################################################################
##Find Null difference between pre/post treatment mean expression##
###################################################################

library(data.table)
library(tidyverse)
library(ggforce)
library(foreach)

#read in allele data for SNP, expression data, and list of genes to be examined
setwd("/home/bal7cg/R/data.objects/")
Alleles <- readRDS("allelesecondary")
#this is the low expression-filtered out data. 
expression_matrix_filename <- readRDS("highcountmatrix")
Genelist <- expression_matrix_filename[,1]
genes414 <- readRDS("genes414")

colnames(genes414)[1] <- "GeneID"
shortplot <- genes414[c(1:3),]

####################################################################
##Run a permutation loop to find read means of different gene sets##
####################################################################
#create vector to determine number of permutations
f = c(1:1000)
#set foreach loop
out <- foreach(f=f)%do%{
#randomly pull out 414 genes
generandom <-  Genelist[sample(nrow(Genelist), 414), ]
#filter out only the 414 genes chosen, created long form data set
Genemelt <- merge(genes414, expression_matrix_filename, by = "GeneID")
Genemelt= melt(Genemelt,id.vars = c("GeneID"),variable.name= "sampleids", value.name = "reads")
Genemelt = Genemelt[, c("pre_pos","cage","line") := tstrsplit(sampleids,"[.]")]
#Turn pre/post names into 0/1 factor values
Genemelt$pre_pos <- str_replace_all(Genemelt$pre_pos,"4CB","0")
Genemelt$pre_pos <- str_replace_all(Genemelt$pre_pos,"23HR","1")
Genemelt$pre_pos <- as.factor(Genemelt$pre_pos)
#find average reads for pre and pos groups
pre.reads <- Genemelt %>% 
  filter(pre_pos == 0)
post.reads <- Genemelt[pre_pos == 1]
pre.summary <- summary(pre.reads$reads)
means <- data.frame(Permutation = f,
                    pre.variance = var(pre.reads$reads), 
                    pre.means = mean(pre.reads$reads),
                    post.variance = var(post.reads$reads), 
                    post.means = mean(post.reads$reads),
                    mean.difference = (mean(post.reads$reads)- mean(pre.reads$reads))
                    )
}
mean.effect.permutations <- rbindlist(out)
#######################################
##Use Bootstrap hypothesis testing##
####################################
#data used: x = whether a given gene is part of 414 group(SNP) or not (NOT)
# y = difference between mean pre starvation gene expression and mean post starvation gene expression
#null hypothesis: whether any given gene is part of the group or not will have no impact on the y value
#alterate hypothesis: genes in the 414 group will have statistically lower y values then others
#test statistic  abs(mean(Y$SNP) - mean(Y$NOT))

########################################
##create the y values described above##
######################################
#recreate the long form data set
Genemelt= melt(expression_matrix_filename,id.vars = c("GeneID"),variable.name= "sampleids", value.name = "reads")
Genemelt = Genemelt[, c("pre_pos","cage","line") := tstrsplit(sampleids,"[.]")]
#Turn pre/post names into 0/1 factor values
Genemelt$pre_pos <- str_replace_all(Genemelt$pre_pos,"4CB","0")
Genemelt$pre_pos <- str_replace_all(Genemelt$pre_pos,"23HR","1")
Genemelt$pre_pos <- as.factor(Genemelt$pre_pos)
#find average reads for pre and pos groups
pre.reads <- Genemelt %>% 
  filter(pre_pos == 0) %>% 
  group_by(GeneID) %>% 
  summarise(means.pre = mean(reads))

post.reads <- Genemelt %>% 
  filter(pre_pos == 1) %>% 
  group_by(GeneID) %>% 
  summarise(means.post = mean(reads))
Gene.means <- merge(pre.reads, post.reads, by = "GeneID")
Gene.means <- Gene.means %>% 
  mutate( mean.dif = means.pre - means.post)
######################################
##create the X values described above#
######################################

Group.decide <- match(Gene.means[[1]], genes414[[1]])
Group.decide <- ifelse(is.na(Group.decide) == TRUE, "NOT", "SNP")
Gene.means <- data.table(
  Gene.means,
  Group.decide
)
head(Gene.means)
#Gene.means column 4 is our Y value, column 5 is the X value. 
###############################################
##Make test statistics for hypothesis testing##
###############################################
test.stat = abs(mean(Gene.means$mean.dif[Group.decide == "SNP"]- mean(Gene.means$mean.dif[Group.decide == "NOT"])))
#now our goal is to create hundreds of bootstrap test statistics with genes randomly assigned to the two groups
#This is under the asssumption that H0 is true, and the groups have no significance.
######################
##set up bootstrap####
######################
#set.seed(112358)# commands to take the same random sample, for reproducibility
set.seed(112)
#set n, the size of our sample 
n  <- length(Gene.means$Group.decide)
#set B, the number of bootstrap resamples
B <- 10000
Y <- Gene.means$mean.dif # the variable that will be sampled
#Create the Boostrap data!
Bootstrapsamples <- matrix(sample(Y, size = n*B, replace = T), nrow = n, ncol = B)
#create an empty vector to store the test statistics
test.vect <- rep(0, B)
#use a For loop to populate this vector with real test stats
for(i in 1:B) {
  #calculate the test stat
  test.vect[i] <- abs(mean(Bootstrapsamples[1:414,i]) #notice that we assign the first 414 elements to the "SNP" group
                      -mean(Bootstrapsamples[415:9361,i])) # and the rest to the "NOT" group
}
head(test.vect, 50
     )
################################
##find p value from test stats##
################################
#remember that the p value is the chance that our observed test statistic could be true given the null hypothesis
#in Bootstrapping, P is determined by the #of bootstrap test stats >= observed test stat, divided by B

(test.vect >= test.stat)[1:20]
#as an exampe, the first 20 boot strap test stats are all lower

mean(test.vect >= test.stat)
#save results
saveRDS(Bootstrapsamples,"Bootstrapsamples")
#what does a simple t test look like here?
t.test(Gene.means$mean.dif[Group.decide == "NOT"], Gene.means$mean.dif[Group.decide == "SNP"])
