#setting the working directory
setwd("F:\\r_projects\\thesis\\1000G_snpr\\heterozygosity\\LDdata")
#vector contains all the family/population names 
i <- c("GBR", "ITU", "FIN", "BEB", "CLM", "CDX", "STU", "MSL", "PJL", "ESN", "ACB", "PUR")

#loading dependencies
library(dplyr)
library(stringr)
library(ggplot2)

# this for loop produces LD values for every 1000kb of genome length and stores in dataframe object 
for (x in i){
    
    dfr <- read.delim(str_interp("LD_1kg_${x}.ld"),sep="",header=T,stringsAsFactors=F)
    
    data1 <- c(dfr$BP_B - dfr$BP_A)
    
    data2 <- data.frame(dfr$R2)
    
    data3 <- data.frame(data1, data2)
    
    data3$distc <- cut(data3$data1,breaks=seq(from=min(data3$data1)-1,to=max(data3$data1)+1,by=1000))
    
    
    dfr1 <- data3 %>% group_by(distc) %>% summarise(mean=mean(dfr.R2),median=median(dfr.R2))
    
    
    dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                            end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                            mid=start+((end-start)/2))
    
    
    dfr1 <- dfr1[-nrow(dfr1),]
    assign(paste0("p_",x), dfr1)
    
}


# this initial code is for the plotting of the LD for different population. This part of the script is primitive. I will fix this later.
ggplot()+

     geom_smooth(data=p_ACB,aes(x=start,y=mean,,colour="grey40"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_BEB,aes(x=start,y=mean,colour="red4"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_CDX,aes(x=start,y=mean,colour="slateblue4"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_CLM,aes(x=start,y=mean,colour="mediumblue"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_ESN,aes(x=start,y=mean,colour="yellow"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_FIN,aes(x=start,y=mean,colour="seagreen2"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_ITU,aes(x=start,y=mean,colour="sienna"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_MSL,aes(x=start,y=mean,colour="green"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_PJL,aes(x=start,y=mean,colour="indianred1"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_PUR,aes(x=start,y=mean,colour="darkgreen"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_STU,aes(x=start,y=mean,colour="cyan"), method = "loess",  size=1.0,alpha=0.5)+
     geom_smooth(data=p_GBR,aes(x=start,y=mean,colour="darkgoldenrod1"), method = "loess",  size=1.0,alpha=0.5)+

    
    labs(x=expression(Distance~(10^{5})),y=expression(LD~(r^{2})))+
    scale_x_continuous(breaks=c(0,2*10^5,4*10^5,6*10^5,8*10^5),labels=c("0","2","4","6","8"))+
    
    scale_color_identity(name = "Population",
                         breaks = c("grey40", "red4", "slateblue4","mediumblue", "yellow", "seagreen2","sienna", "green", "indianred1","darkgreen", "cyan", "darkgoldenrod1"),
                         labels = c("ACB", "BEB", "CDX", "CLM", "ESN","FIN", "ITU", "MSL", "PJL", "PUR", "STU", "GBR"),
                         guide = "legend") +
    theme_bw()
