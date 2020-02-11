# Manon Morin, November 30th 2018
#That script compares the gene fitness values calculated by wetmore approach amd our approach

rm(list=ls())

  #### Packages and functions ####

library(tidyr)
library(dplyr)
library(ggplot2)
library(gridExtra)

  #### Data upload ####

load("Averaged_Replicates_HOI_Ecoli.Rdata")
Data_new=Final_gene_Fitness

Data_Wet=read.csv("Fitness_Wetmore.csv")

  #### Filtering of genes with fitness in both approaches ####
Genes_us=unique(Data_new$sysName)

Genes_Names=intersect(Data_Wet$sysName, Genes_us)

Data_new1=subset(Data_new, Data_new$sysName%in%Genes_Names)
Data_Wet1=subset(Data_Wet, Data_Wet$sysName%in%Genes_Names)

  ### Distribution of fitness, Wetmore + PCA ###

  # Distribution
Data_Wet_dis=Data_Wet1%>%gather(Cdt, Fit, Alone:Com)
plot_Wet=ggplot(Data_Wet_dis, aes(x=Fit, col=Cdt)) + geom_density(aes(fill=Cdt)) + theme_light() + ggtitle("Distribution of Wetmore fitness") +
  facet_grid(Cdt ~ .) + labs(y = "Density", x="Wetmore fitness")

plot_new1=ggplot(Data_new1, aes(x=CenteredFit, col=Cdt)) + geom_density(aes(fill=Cdt)) + theme_light() + ggtitle("Distribution of Wetmore fitness") +
  facet_grid(Cdt ~ .) + labs(y = "Density", x="Wetmore fitness")

 # PCA 1 - Eigen

Data_PCA=Data_Wet1[,2:9]
row.names(Data_PCA)=Data_Wet1$sysName

PCA_eigen=prcomp(Data_PCA, scale=TRUE)
plot_PCA1=fviz_pca_var(PCA_eigen,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE )    # Avoid text overlapping

 # PCA 2 - Points

Data_plot2=t(Data_PCA)

Data_PCA2=Data_plot2

Data_plot2=as.data.frame(Data_plot2)
Data_plot2$Cdt=rownames(Data_plot2)

plot_PCA2=autoplot(prcomp(Data_PCA2), data = Data_plot2, colour="Cdt", size=5)+theme_light()

  #### Correlation between both approaches ####

Data_both=data.frame()

Cdt_list=c(unique(as.character(Data_new1$Cdt)))

for (i in 1:length(Cdt_list)){
  Dat1=subset(Data_new1, Data_new1$Cdt==Cdt_list[i])
  Dat2=subset(Data_Wet_dis, Data_Wet_dis$Cdt==Cdt_list[i])
  DatF=left_join(Dat1, Dat2, by="sysName")
  
  Data_both=rbind(Data_both, DatF)
}


 # Obtention of correlation score
library(DescTools)
Table_CorP=data.frame("Cdt"=Cdt_list,"Cor"=0)

for (i in 1:length(Cdt_list)){
  Dat=subset(Data_both, Data_both$Cdt.x==Cdt_list[i])
  corV=cor(Dat$CenteredFit, Dat$Fit)
  Table_CorP[i,2]=corV
}

Table_CorS=data.frame("Cdt"=Cdt_list,"Cor"=0)

for (i in 1:length(Cdt_list)){
  Dat=subset(Data_both, Data_both$Cdt.x==Cdt_list[i])
  corV=cor(Dat$CenteredFit, Dat$Fit, method=c("spearman"))
  Table_CorS[i,2]=corV
}

Table_CorC=data.frame("Cdt"=Cdt_list,"Cor"=0)

for (i in 1:length(Cdt_list)){
  Dat=subset(Data_both, Data_both$Cdt.x==Cdt_list[i])
  corC=CCC(Dat$CenteredFit, Dat$Fit)
  Table_CorC[i,2]=corC$rho.c$est
}

List_plot=list()

for (i in 1:length(Cdt_list))
  local({
  Dat=subset(Data_both, Data_both$Cdt.x==Cdt_list[i])
  plot=ggplot(Dat, aes(Dat$CenteredFit, Dat$Fit)) + geom_point() + theme_light() +
    ggtitle(paste(Cdt_list[i])) + labs(x = "New fitness", y="Wetmore fitness") + stat_smooth(method="lm", se=FALSE, color="blue") + 
    geom_abline(intercept = 0, slope = 1, color="red", linetype="dashed", size=1.5)
  
  List_plot[[i]]<<-plot
})



grid.arrange(List_plot[[1]],List_plot[[2]],List_plot[[3]],List_plot[[4]],
             List_plot[[5]],List_plot[[6]],List_plot[[7]],List_plot[[8]], nrow=2, ncol=4)

grid.arrange(List_plot[[1]]:List_plot[[8]], nrow=2, ncol=4)












