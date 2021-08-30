#setting working directory
setwd("F:\\r_projects\\thesis\\1000G_snpr\\heterozygosity")

#vector contains the population names
i <- c("GBR", "ITU", "FIN", "BEB", "CLM", "CDX", "STU", "MSL", "PJL", "ESN", "ACB", "PUR")

# produces vector for every population containing number of alleles in the interval of .1 from 0 to .5
for (x in i) {
    s <- str_interp("plink2 --bfile ${x}_1kg --freq --out allele_frequency_${x}")
    system(s)
    
    any_freq <- assign(paste0(x,"_freq"), read.csv(paste0("allele_frequency_",x,".afreq"), sep = '\t', header = TRUE))
    any_vec <- c()

    h=1
    while (h != 6) {
        
        gen_var_x <- nrow(any_freq %>% filter(ALT_FREQS <= h/10, ALT_FREQS > ((h/10)-.1)))
        any_vec[h] <- gen_var_x
        h <- h + 1
        
    }
    assign(paste0(x,"_vec"), any_vec)
    rm(gen_var_x)
    rm(any_freq)
    rm(any_vec)
    print("yeas")
}



### making all the population vector object stratified 

tr<- c()

for (x in i){
    vec_a <- get(str_interp("${x}_vec"))
    tr <- rbind(tr, vec_a)
}

rownames(tr) <- i
colnames(tr) <- c("tr1",'tr2','tr3', 'tr4', 'tr5')



### default barplot from r

barplot(tr$tr5, names.arg = i, xlab = "population", las = 1, ylab = "number of allele", ylim = c(0,40000), main = "0.4< Allele Frequency >=0.5",
        col = brewer.pal(12, name = "Paired"), border  = brewer.pal(12, name = "Paired"))



### plotting with ggplot

ggplot(data = tr) +
    geom_bar(mapping = aes(y = tr1, x = rownames(tr)), stat = "identity",
             col = brewer.pal(12, name = "Paired"), fill = brewer.pal(12, name = "Paired"))+
    labs(x = "number of populations", y = "allele count")

     