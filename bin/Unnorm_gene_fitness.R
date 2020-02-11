##### UN-NORMALIZED GENE FITNESS CALCULATION #####

#This function first calculates strain fitness and then un-normalizedgene fitness

Unnorm_gene_fitness=function(Data,cdnbF,genes.tab,org_locId){

# First we calculate the log2 value of each count

  Data_log=Data

  for (i in 8:dim(Data_log)[2]){
    Data_log[,i]=log2(Data_log[,i])
  }

  #### (i) Calculating insertion mutant (strain) fitness ####
  
#To get the strain fitness we do log2condition - logT0
#We first isolate the T0 columns that we will use for all conditions
#Then we substract the corresponding T0 values to each condition - Before substracting we want to make surecthat both tables are sorted the same way, so we order the dat based on the pos column 

  T0_val=Data_log[,-9:-(ncol(Data_log))]
  Cdt_val=Data_log[,-8]

  T0_val=T0_val[order(T0_val$pos,T0_val$barcode),]
  Cdt_val=Cdt_val[order(Cdt_val$pos, Cdt_val$barcode),]

  StrainFit_Table=Cdt_val

  for (i in 8:dim(StrainFit_Table)[2]){
    StrainFit_Table[,i]=(Cdt_val[,i]-T0_val[,8])
  }
  
  #### (ii) Calculating raw gene fitness #####
  
  # The raw gene fitness is the average of all the strain fitness in that gene
  
  if (org_locId=="Num"){
    Locus=unique(as.numeric(StrainFit_Table$locusId))
  }
  
  if (org_locId=="Char"){
    Locus=unique(as.character(StrainFit_Table$locusId))
  }
  
  Mat_GeneFit_mean=matrix(0,nrow=length(Locus), ncol=(cdnbF))
  colnames(Mat_GeneFit_mean)=colnames(StrainFit_Table[,8:ncol(StrainFit_Table)])
  rownames(Mat_GeneFit_mean)=Locus
  
  for (i in 8: dim(StrainFit_Table)[2]){
    Data_Cdt=data.frame(StrainFit_Table[,1:7],StrainFit_Table[i])
    for (j in 1:length(Locus)){
      Loc=Locus[j]
      Temp=subset(Data_Cdt, Data_Cdt$locusId==Loc)
      Mat_GeneFit_mean[j,(i-7)]=mean(Temp[,8])
    }
  }
  
  Table_GeneFitness_MEAN=data.frame("locusId"=Locus,Mat_GeneFit_mean)  
  
  #We do a similar procedure to obtain the variance of gene fitness
  
  
  Mat_GeneFit_var=matrix(0,nrow=length(Locus), ncol=(cdnbF))
  colnames(Mat_GeneFit_var)=colnames(StrainFit_Table[,8:ncol(StrainFit_Table)])
  rownames(Mat_GeneFit_var)=Locus
  
  for (i in 8: dim(StrainFit_Table)[2]){
    Data_Cdt=data.frame(StrainFit_Table[,1:7],StrainFit_Table[i])
    for (j in 1:length(Locus)){
      Loc=Locus[j]
      Temp=subset(Data_Cdt, Data_Cdt$locusId==Loc)
      Mat_GeneFit_var[j,(i-7)]=var(Temp[,8])
    }
  }
  
  Table_GeneFitness_VAR=data.frame("locusId"=Locus,Mat_GeneFit_var)
 
#We want to add the sysName and scaffold ID and gene position
  
  Table_GeneFitness_MEAN=Table_GeneFitness_MEAN[order(Table_GeneFitness_MEAN$locusId),]
  Table_GeneFitness_VAR=Table_GeneFitness_VAR[order(Table_GeneFitness_VAR$locusId),]
  
  genes.tab=genes.tab[order(genes.tab$locusId),]
  genes.tab=subset(genes.tab, genes.tab$locusId%in%Table_GeneFitness_MEAN$locusId)
  
  
  Table_GeneFitness_MEAN$sysName=genes.tab$sysName
  Table_GeneFitness_MEAN$begin=genes.tab$begin
  
  Tables_unnorm_list=list("Strain_fitness"=StrainFit_Table,"Gene_unnorm_fitness"=Table_GeneFitness_MEAN,
                  "Gene_fitness_variance"=Table_GeneFitness_VAR)
  Tables_unnorm_list
}