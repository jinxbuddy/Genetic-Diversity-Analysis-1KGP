#setwd

setwd("F:\\r_projects\\thesis\\1000G_snpr")

#load dependencies
library(ggplot2)
library(LEA)
library(mapplots)
library(fields)
library(tidyverse)
# selection of a subset of data by families

read_tsv("1000G_phase3_common_norel.fam", col_names = F) %>%
    select(X1,X2) %>%
    filter(X1 == "GBR" | X1 == "FIN" | X1 == "PUR" | X1 == "PJL" | X1 == "CDX" | X1 == "ACB" | X1 == "ESN" | X1 == "BEB" | X1 == "STU" | X1 == "ITU" | X1 == "MSL" | X1 == "CLM") %>%
    write_delim("1kg_hapmap3_12pop.txt", col_names = FALSE)



#

system("plink --bfile 1000G_phase3_common_norel --keep 1kg_hapmap3_12pop.txt --make-bed --out 12pop_1kg_hapmap3")

#missingness filters

system("plink --bfile 12pop_1kg_hapmap3 --geno 0.1 --mind 0.25 --threads 2 --make-bed --out missingness_filtered_data ")

#minor allele frequency reporting

system("plink --bfile missingness_filtered_data --freq --out alle_frequency_plink1.9")

system("plink --bfile missingness_filtered_data --threads 2 --freqx --out allele_frequency_plink1.9")

system("plink2 --bfile missingness_filtered_data --freq alt1bins=0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4 --threads 2 --out allele_spectrum")


#Selectinga Sample Subset Without Very Close Relatives

system("plink2 --bfile missingness_filtered_data --king-cutoff 0.177 --threads 2 --make-bed --out relpruned_data")




#hardy-weinberg


system("plink2 --bfile missingness_filtered_data --hwe 1e-25 keep-fewhet --make-bed --out hwe_filtered_data")


# Selecting a SNP Subset in Approximate Linkage Equilibrium


system("plink --bfile hwe_filtered_data --indep-pairwise 200 50 0.2 --threads 2 --out ldpruned_snplist")

system("plink --bfile hwe_filtered_data --extract ldpruned_snplist.prune.in --threads 2 --make-bed --out ldpruned_data")


#maf filtered data

system("plink --bfile ldpruned_data --maf 0.1 --threads 2 --make-bed --out maf_filtered_data")

#PCA

system("plink2 --bfile maf_filtered_data --pca --threads 2 --out pca_results")


#loading the data as a list

pca_table <- read.table("pca_results.eigenvec", header=TRUE, comment.char="")

#loading the ggplot library

library(ggplot2)

# PCA plot

ggplot(data = pca_table) +
    geom_point(mapping = aes(x = PC2, y = PC3, color = X.FID, ), size = 2,  show.legend = TRUE ) + 
    geom_hline(yintercept = 0, linetype="dotted") + 
    geom_vline(xintercept = 0, linetype="dotted") +
    
    labs(title = "PCA Plot") + 
    theme_minimal()


##### the admixture analysis


#converting ped to vcf

system("plink --bfile maf_filtered_data --recode vcf --threads 2 --out 12_pop_1kg")

#converting vcf to lfmm
vcf2lfmm("12_pop_1kg.vcf")

# convertinf lfmm to geno
lfmm2geno("12_pop_1kg.lfmm")

# running admixture
project3.snmf = snmf("12_pop_1kg.geno",
                     K = 13:15, 
                     entropy = TRUE, 
                     repetitions = 10, 
                     project = "new")

# plot cross-entropy criterion of all runs of the project
plot(project3.snmf, cex = 2, col = "red3", pch = 19)

# get the cross-entropy of the 10 runs for K = 4
ce3 = cross.entropy(project3.snmf, K =12 )

# select the run with the lowest cross-entropy for K = 4
best2 = which.min(ce3)

# display the Q-matrix

my.colors <- c("red4", "slateblue4", 
               "mediumblue", "yellow", "seagreen2", "sienna","green", "indianred1", "dimgray", "darkgreen","cyan","darkgoldenrod1")


barchart(project3.snmf, K = 12, run = best2, 
         border = NA, space = 0, col = my.colors, 
         xlab = "Individuals", ylab = "Ancestry proportions", 
         main = "Ancestry matrix") -> bp

axis(1, at = 1:length(bp$order), 
     labels = bp$order, las = 3, cex.axis = .4)





###For further data exploration

# show the project
show(project2.snmf)

# summary of the project
summary(project2.snmf)

# get the cross-entropy for all runs for K = 4
ce7 = cross.entropy(project2.snmf, K = 7)

# get the cross-entropy for the 2nd run for K = 4
ce4 = cross.entropy(project2.snmf, K = 4, run = 2)

# get the ancestral genotype frequency matrix, G, for the 2nd run for K = 4. 
freq = G(project2.snmf, K = 7, run = 2)


# display the Q-matrix
Q.matrix <- as.qmatrix(Q(project2.snmf, K = 7, run = best2))
my.colors <- c("red4", "slateblue4", 
               "mediumblue", "yellow", "seagreen2", "sienna","green", "indianred1", "dimgray", "darkgreen")

barplot(Q.matrix, 
        border = NA, 
        space = 0, 
        col = my.colors, 
        xlab = "Individuals",
        ylab = "Ancestry proportions", 
        main = "Ancestry matrix") -> bp

axis(1, at = 1:nrow(Q.matrix), labels = bp$order, las = 3, cex.axis = .4)



