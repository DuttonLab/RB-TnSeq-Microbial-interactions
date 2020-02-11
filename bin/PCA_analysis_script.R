## Manon Morin- November 30th

## That script is part of the RB-TnSeq analysis to generate PCA plots of the conditions
rm(list=ls())

#Packages and functions

library(ggplot2)
library(factoextra)
library(tidyr)
library(dplyr)

    #### Data upload ####

# Data upload
#You can upload the Rdat file saved at the end of the PartII : Averaging_replicates.R

load("Averaged_Replicates_HOI_Ecoli.Rdata")

DataAll=Final_gene_Fitness   # For PCA on average

Data_byRep=AllReplicate  #For PCA using replicates


    #### I) PCA using fitness average values #####

#Here we will do 2 different PCA: i) plotting the eigen vectors to see how each condition drives gene fitness values and ii) the PCA to see on overall how the conditions are diverging/converging based on global gene fitness

# 1) Plotting eigen vectors

  # To plot the eigen vector we need a table where each column is a condition and each row a gene

Data_plot=data.frame("sysName"=DataAll$sysName,"Fit"=DataAll$WeightedVar,"Cdt"=DataAll$Cdt)
Data_plot=Data_plot %>% spread(Cdt, Fit)

Data_PCA=Data_plot[,2:9]
row.names(Data_PCA)=Data_plot$sysName

PCA_eigen=prcomp(Data_PCA, scale=TRUE)
fviz_eig(PCA_eigen)
plot_eigen=fviz_pca_var(PCA_eigen,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# 2) Plot the condition

Data_plot2=t(Data_PCA)

Data_PCA2=Data_plot2

Data_plot2=as.data.frame(Data_plot2)
Data_plot2$Cdt=rownames(Data_plot2)

PCA_cdt=prcomp(Data_PCA2, scale=TRUE)
fviz_eig(PCA_cdt)

plot_PCA=autoplot(prcomp(Data_PCA2), data = Data_plot2, colour="Cdt", size=5)+theme_light()


    #### I) PCA using fitness replicates values #####

# We do the same PCA than before but including replicates
Data_plotR=data.frame("sysName"=Data_byRep$sysName,"Fit"=Data_byRep$Fitness_byMean,"Cdt"=Data_byRep$Cdt, "Rep"=Data_byRep$Rep)
Data_plotR= Data_plotR %>% unite(CdtR, Cdt,Rep, sep = ".")
Data_plotR=Data_plotR %>% spread(CdtR, Fit)

Data_PCA=Data_plotR[,2:25]
row.names(Data_PCA)=Data_plotR$sysName

PCA_eigen=prcomp(Data_PCA, scale=TRUE)
#fviz_eig(PCA_eigen)
plot_eigenRep=fviz_pca_var(PCA_eigen,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

# 2) Plot the condition

Data_plot3=t(Data_PCA)

Data_PCA3=Data_plot3

Data_plot3=as.data.frame(Data_plot3)
Data_plot3$Cdt=c("Alone","Alone","Alone","Com","Com","Com", "Geo","Geo","Geo","GeoSam","GeoSam","GeoSam","Haf","Haf","Haf","HafGeo",
                 "HafGeo", "HafGeo","HafSam","HafSam","HafSam","Sam","Sam","Sam")
Data_plot3$Rep=c("R1","R2","R3")

plot_PCARep=autoplot(prcomp(Data_PCA3), data = Data_plot3, colour="Cdt", size=5, shape="Rep")+theme_light()





