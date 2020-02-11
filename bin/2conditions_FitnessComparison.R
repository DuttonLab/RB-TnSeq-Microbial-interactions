############################  GENE FITNESS COMPARISON BETWEEN TWO CONDITIONS  ############################
############################         PART III: RB-TNSEQ comparative analysis        ############################

# Compares genes fitness values of all conditions against a chosen reference condition
# Append a "Category" to each compared gene mentionning if fitness values are significantly different or not

rm(list=ls())

# Package and functions upload

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)

source("Category_definition.R")
source("Comparison_test.R")

# Step 1: Data import and parameters settings
  # Import the .Rdata file generated in Averaging_replicates.R
  # Set up org_locId
  # Set up your condition of reference (has to be one of your conditions)
  # Set up your alpha value for the T-test

load("Averaged_Replicates_HOI_Ecoli.Rdata")   # here you upload the .RData ouput of Averaging_replicates.R
Data_Fitness=Final_gene_Fitness

org_locId="Num"  # Again, you set up your organisms whether "Ecoli" if the locusId are in numeric form or "Pseudo" if the locusId are in character form
Condition1="Alone"  #here you write the name of one of the 2 conditions you want to compare. Make sure it is the same name than previously used
alpha=0.002   #here you choose any alpha you want. Just be aware that it is where you can control the amount of false discovery you allow

# Step 2: Run the comparison for all conditions against the reference one
  # All the conditions that are not the reference are going to be tested against the refence one after the other using a for loop embedding the comparison function


Table_all=data.frame() # creates a table to eventually store the data
List_plots=list(0,0,0,0,0,0,0) # creates a list to store comparison of plots for each comparison in the "for loop" . 
                               # needs to have as many spots in the list than you have comparison

Conditions=as.vector(unique(Data_Fitness$Cdt))
Conditions=Conditions[-1] # generates a vector with all the conditions to compare against the reference one


for (i in 1:length(Conditions)){
  Cdt=Conditions[i]
  Condition2=paste(Cdt)  
  Test1=Comparison_test(Data_Fitness,Condition1, Condition2, alpha, org_locId)   
  #You obtain a list
  #list[[1]]= the table containing all the genes tested + fitness + variance values in both conditions + Fcalc (Fisher value calculates for the test of variance + Tcalc (Tscore value for the student test))
  # and the category column indicates if the difference of fitness is significant (Sig), if it is not (Not_Sig) or if the student was not performed because of unequal variances (Not_tested)
  #list[[2]]=the scatter plot representing the data
  
  Table_all=rbind(Table_all, Test1[[1]])
  List_plots[[i]]= Test1[[2]]
}

#Visualize all the scatter plots
grid.arrange(List_plots[[1]],List_plots[[2]],List_plots[[3]],List_plots[[4]],
             List_plots[[5]],List_plots[[6]],List_plots[[7]], nrow=3, ncol=3)

write.csv(Table_all, "")
save(Table_all, List_plots, file=".RData")

