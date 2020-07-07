##### DATA FORMATING IN CASE OF CONDITIONS USING DIFFERENT T0s #####
### Replacement for Step1 of Averaging_replicates2.R

#If you are comparing experiments with different T0's, you will need to combine the All_data_Replicate tables from the separate tables from Script I prior to binding replicates together. 
#This can easily be done using code similar to this for each replicate

#This code is written if you have 2 sets of conditions with different 2 different T0s
#For each replicate you'll have to import the .Rdata file of each set of conditions (generated in the Gene_Fitness_Replicate.R)
#For each replicate, you need to:
  #-combine the fitness values as a single table
  #-combine the fitness mean values as a single vector

##Replicate1:
load("") 
Replicate1=All_data_Replicate

load("") 
Replicate1x=All_data_Replicate

locus_shared=intersect(Replicate1$locusId,Replicate1x$locusId) # Because using different T0 can lead to different loci set, we have to select the loci that are shared in all conditions
Replicate1=subset(Replicate1, Replicate1$locusId%in%locus_shared)
Replicate1x=subset(Replicate1x, Replicate1x$locusId%in%locus_shared)

Replicate1=rbind(Replicate1,Replicate1x)
head(Replicate1,n=5)

##Replicate2:
load("") 
Replicate2=All_data_Replicate


load("") 
Replicate2x=All_data_Replicate

locus_shared=intersect(Replicate2$locusId,Replicate2x$locusId)
Replicate2=subset(Replicate2, Replicate2$locusId%in%locus_shared)
Replicate2x=subset(Replicate2x, Replicate2x$locusId%in%locus_shared)

Replicate2=rbind(Replicate2,Replicate2x)
head(Replicate2,n=5)

##Replicate3:
load("") 
Replicate3=All_data_Replicate

load("") 
Replicate3x=All_data_Replicate

locus_shared=intersect(Replicate3$locusId,Replicate3x$locusId)
Replicate3=subset(Replicate3, Replicate3$locusId%in%locus_shared)
Replicate3x=subset(Replicate3x, Replicate3x$locusId%in%locus_shared)

Replicate3=rbind(Replicate3,Replicate3x)
head(Replicate3,n=5)

# We bind all replicates together and add a column "Rep" to identify were it is coming from

Replicate1$Rep="R1"
Replicate2$Rep="R2"
Replicate3$Rep="R3"

AllReplicate=rbind(Replicate1,Replicate2,Replicate3)
head(AllReplicate, n=5)

# Save

save(AllReplicate, file="") # save your table and mean vectors to be imported in step 2 of "Averaging_replicate2.R"
