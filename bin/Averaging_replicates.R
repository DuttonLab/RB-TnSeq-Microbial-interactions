############################  GENE FITNESS CALCULATION: AVERAGING ACROSS REPLICATES  ############################
############################         PART II: RB-TNSEQ comparative analysis        ############################

# Averages gene fitness across replicates
# Requires .RData files generated for each replicate in the script "Gene_Fitness_Replicate.R"
# Produces a final file file containing final normalized gene fitness and associated fitness variance across replicates 
# Many plots are generated along data processing to visually follow data transformation

# Step 0: Packages and source functions
rm(list=ls())

library(dplyr)
library(ggplot2)
library(ggrepel)
library(DescTools)

source("Correlation_Rep.R")
source("Weighted_average.R")
source("Mean_Correction.R")
source("Mode_Centered.R")

# Step1: Data import and parameters set up
 # Import .Rdata files generated for each replicate in Gene_Fitness_Replicate.R
 # Isolate and rename All_data_Replicate just after loading to avoid overwritting.

load("HOI_Ecoli1.Rdata") 
Replicate1=All_data_Replicate
head(Replicate1,n=5)
mean1=Data_norm_mean_OP[[2]]

load("HOI_Ecoli2.Rdata") 
Replicate2=All_data_Replicate
head(Replicate2,n=5)
mean2=Data_norm_mean_OP[[2]]

load("HOI_Ecoli3.Rdata") 
Replicate3=All_data_Replicate
head(Replicate3,n=5)
mean3=Data_norm_mean_OP[[2]]

  # We bind all replicates together and add a column "Rep" to identify were it is coming from

Replicate1$Rep="R1"
Replicate2$Rep="R2"
Replicate3$Rep="R3"

AllReplicate=rbind(Replicate1,Replicate2,Replicate3)
head(AllReplicate, n=5)

# Set up parameter org_locId
org_locId="Num"   # "Num" if numeric , "Char" if characters

# Step 2: Replicate correlation calculation and visualization (for gene fitness values)
  # Determines the correlation between each replicate, stores them in a matrix 
  # Visualizes correlation in two different ways: (i)  usual correlation plots and (ii) distribution of correlation across all conditions and replicates

  # Plots for fitness values or variance distribution in each condition and replicate
plot1_fit=ggplot(AllReplicate, aes(x=Fitness_byMean,col=Cdt))+geom_density(aes(fill=Cdt)) + theme_light()+
  ggtitle("Distribution of gene fitness values") + facet_grid(Rep ~ Cdt) + labs(x = "Gene fitness value", y="Density")

plot1_var=ggplot(AllReplicate, aes(x=Fitness_Variance,col=Cdt))+geom_density(aes(fill=Cdt)) + theme_light()+
  ggtitle("Distribution of variance") + facet_grid(Rep ~ Cdt) + labs(x = "Variance associated with fitness values", y="Density")


  # Calculation of correlation for all pairs of replicate (and each condition) + plots
Correlation_table_fit=Correlation_Rep(AllReplicate) # Calculates different correlation coefficient + save plots

head(Correlation_table_fit,n=5)

plot_allPcor_fit=ggplot(Correlation_table_fit, aes(x=Cond,y=P_Rsquared)) + geom_point(size=2, aes(col=Comp))+
  theme_light() + labs(y = "Pearson Correlation coefficient", x="Condition")

plot_allScor_fit=ggplot(Correlation_table_fit, aes(x=Cond,y=P_Rsquared)) + geom_point(size=2, aes(col=Comp))+
  theme_light() + labs(y = "Spearman Correlation coefficient", x="Condition")


plot_allLcor_fit=ggplot(Correlation_table_fit, aes(x=Cond,y=P_Rsquared)) + geom_point(size=2, aes(col=Comp))+
  theme_light() + labs(y = "Lin's Correlation coefficient", x="Condition")


# Step 3: Weighted average of gene fitness across replicate
 # Averages gene fitness values across replicates for each condition
 # Visualizes fitness values or variance ditribution per condition 

Average_fitness=Weighted_average(AllReplicate, org_locId)

plot2_fit=ggplot(Average_fitness, aes(x=WeightedFit,col=Cdt))+geom_density(aes(fill=Cdt)) + theme_light()+
  ggtitle("Distribution of average fitness values across replicates") + facet_grid(Cdt ~ .) + labs(x = "Gene fitness values", y="Density")


plot2_var=ggplot(Average_fitness, aes(x=WeightedVar,col=Cdt))+geom_density(aes(fill=Cdt)) + theme_light()+
  ggtitle("Distribution of variance") + facet_grid(Cdt ~ .) + labs(x = "Variance associated with fitness values", y="Density")


# Step 4: Correction of fitness values by the Mean of means
  #Normalizes each gene fitness value by the Mean of the replicates distribution means (mean1, mean2 and mean3)
  #Visualize fitness values per condition 
Mean_corrected_averageF=Mean_Correction(Average_fitness, mean1,mean2,mean3, org_locId)

plot3_fit=ggplot(Mean_corrected_averageF, aes(x=CorrectedFit,col=Cdt))+geom_density(aes(fill=Cdt)) + theme_light()+
  ggtitle("Distribution of average fitness values across replicates (Mean corrected)") + facet_grid(Cdt ~ .) + labs(x = "Gene fitness values", y="Density")


# Step 5: Center distribution around the mode so neutral fitness=0
 #Normalizes each gene fitness value by the mode of the distribution
 #Visualize fitness values per condition 

Mean_centered_average=Mode_Centered(Mean_corrected_averageF,org_locId)

plot4_fit=ggplot(Mean_centered_average, aes(x=CenteredFit,col=Cdt))+geom_density(aes(fill=Cdt)) + theme_light()+
  ggtitle("Distribution of average fitness values across replicates (Corrected and centered)") + facet_grid(Cdt ~ .) + labs(x = "Gene fitness values", y="Density")


# Step6: Saving data
if(org_locId=="Num"){
  Final_gene_Fitness=left_join(genes.tab,Mean_centered_average, by=c("locusId")) %>% select(-c(type,strand,GC,nTA))
  Final_gene_Fitness=na.omit(Final_gene_Fitness)
}

if(org_locId=="Char"){
  Mean_centered_average$locusId=as.character(Mean_centered_average$locusId)  # the locusId in the Mean table are factors, we need to turn then into charcater for the left_join
  Final_gene_Fitness=left_join(genes.tab,Mean_centered_average, by=c("locusId")) %>% select(-c(type,strand,GC,nTA))
  Final_gene_Fitness$name="No_name"   #replace the NA by something else, otherwise the next NAomit removes everything
  Final_gene_Fitness=na.omit(Final_gene_Fitness)
}

write.csv(Final_gene_Fitness,"")

save(Final_gene_Fitness, Mean_corrected_averageF, genes.tab,
     Average_fitness, Correlation_table_fit,AllReplicate,mean1,mean2,mean3,
     plot1_fit,plot1_var,plot_allcor_fit,plot2_fit, plot2_var,plot3_fit,plot4_fit, file="")

