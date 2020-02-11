####### FUNCTION TO COMPARE FITNESS VALUES BETWEEN TWO CONDITIONS  ###########




# STEP1: Test for equal variance
#Student comparison test can only be performed on samples with equal variance, we first have to identify which genes fitness have equal variance between the conditions
#We perform a Fisher test to test if variances are equal with (p-value of 5%)
#Briefly, this involves for each gene, calculating the ratio of variance (higher variance always as the numerator). If that ratio is GREATER or equal to a given Fisher statistics, then variances are not equal and difference for that gene fitness value cannot be evaluated.
#We use the Fisher statistics for alpha=5% (so F0.025 since we do bilateral test) for constant degrees of freedom of 2 and 2 since we have 3 replicates

Comparison_test=function(Data_Fitness, Condition1, Condition2, alpha, org_locId){
  
  Data_1=subset(Data_Fitness, Data_Fitness$Cdt==Condition1)
  Data_2=subset(Data_Fitness, Data_Fitness$Cdt==Condition2)
  
  Data_1=Data_1[order(Data_1$locusId),]
  Data_2=Data_2[order(Data_2$locusId),]
  
  Locus_list=c(unique(Data_1$locusId))
  
  Mat_Fcalc=matrix(0,nrow=length(Locus_list), ncol=1)
  rownames(Mat_Fcalc)=Locus_list
  colnames(Mat_Fcalc)=c("Fcalc")
  
  #First loop to calculate Fcalc (ratio of variance with maximun as numerator) for each gene
  for (i in 1:length(Locus_list)){
    DatC1=subset(Data_1, Data_1$locusId==Locus_list[i])
    DatC2=subset(Data_2, Data_2$locusId==Locus_list[i])
    V1=max(DatC1$WeightedVar,DatC2$WeightedVar)
    V2=min(DatC1$WeightedVar,DatC2$WeightedVar)
    
    Fcalc=V1/V2
    
    Mat_Fcalc[i,1]=Fcalc
    
  }
  
  if (org_locId=="Num"){
    Data_Fcalc=data.frame("locusId"=as.numeric(rownames(Mat_Fcalc)), "Fcalc"=Mat_Fcalc[,1])
  }
  
  if (org_locId=="Char"){
    Data_Fcalc=data.frame("locusId"=as.character(rownames(Mat_Fcalc)), "Fcalc"=Mat_Fcalc[,1])
  }
  
  Data_FcalcGOOD=subset(Data_Fcalc, Data_Fcalc$Fcalc<=39)
  
  # STEP2: Run student test to compare mean associated with equal variance
  # We use 3 replicates 
  # Depending on the alpha chosen, we can use diffrent cutoffs to identify significant diffrences in fitness
  
  Data_1G=subset(Data_1, Data_1$locusId%in%Data_FcalcGOOD$locusId)
  Data_2G=subset(Data_2, Data_2$locusId%in%Data_FcalcGOOD$locusId)
  
  if(org_locId=="Num"){
    Locus_list2=as.numeric(Data_1G$locusId)
  
    Mat_Tcalc=matrix(0,nrow=length(Locus_list2), ncol=1)
    rownames(Mat_Tcalc)=Locus_list2
    colnames(Mat_Tcalc)=c("Tcalc")
  
  
    for (i in 1:length(Locus_list2)){
    
      DatC1=subset(Data_1G, Data_1G$locusId==Locus_list2[i])
      DatC2=subset(Data_2G, Data_2G$locusId==Locus_list2[i])
    
      Var1=as.numeric(DatC1$WeightedVar)
      Var2=as.numeric(DatC2$WeightedVar)
    
      Fit1=as.numeric(DatC1$CenteredFit)
      Fit2=as.numeric(DatC2$CenteredFit)
    
      Sp=sqrt((2*Var1+2*Var2)/4)
      Tcalc=(Fit1-Fit2)/(Sp*sqrt(2/3))
    
      Mat_Tcalc[i,1]=Tcalc
    
    }
    Data_Tcalc=data.frame("locusId"=as.numeric(rownames(Mat_Tcalc)), "Tcalc"=Mat_Tcalc[,1])
  }
  
  if(org_locId=="Char"){
    Locus_list2=as.character(Data_1G$locusId)
    
    Mat_Tcalc=matrix(0,nrow=length(Locus_list2), ncol=1)
    rownames(Mat_Tcalc)=Locus_list2
    colnames(Mat_Tcalc)=c("Tcalc")
    
    
    for (i in 1:length(Locus_list2)){
      
      DatC1=subset(Data_1G, Data_1G$locusId==Locus_list2[i])
      DatC2=subset(Data_2G, Data_2G$locusId==Locus_list2[i])
      
      Var1=as.numeric(DatC1$WeightedVar)
      Var2=as.numeric(DatC2$WeightedVar)
      
      Fit1=as.numeric(DatC1$CenteredFit)
      Fit2=as.numeric(DatC2$CenteredFit)
      
      Sp=sqrt((2*Var1+2*Var2)/4)
      Tcalc=(Fit1-Fit2)/(Sp*sqrt(2/3))
      
      Mat_Tcalc[i,1]=Tcalc
      
    }
    Data_Tcalc=data.frame("locusId"=as.character(rownames(Mat_Tcalc)), "Tcalc"=Mat_Tcalc[,1])
  }
  
  Cutoff=abs(qt(alpha/2,4))
  
  Data_TcalcGOOD=subset(Data_Tcalc, abs(Data_Tcalc$Tcalc)>=Cutoff)
  
  #Now Data formatting
  
  Data_2C=Data_2[,-2:-7]
  
  Data_test=left_join(Data_1,Data_2C, by=c("locusId"))
  colnames(Data_test)[8:13]=c("Fit_C1","Var_C1","Cond1","Fit_C2","Var_C2","Cond2")
  
  Data_test=left_join(Data_test, Data_Fcalc, by=c("locusId"))
  Data_test=left_join(Data_test, Data_Tcalc, by=c("locusId"))
  
  Data_test=Category_definition(Data_test, alpha)  #here you can change the value of alpha you want to use (0.05, 0.01, 0.005, 0.001)
  
  
  scatter_plot=ggplot(Data_test %>% arrange(Category),aes(x=Fit_C1,y=Fit_C2,color=Category))+geom_point(aes(color=Category, size=Category)) +
    geom_abline(intercept = 0, slope = 1, color="red", size=0.5) + scale_color_manual(values=c("Gray70","Gray95","blue"))+
    scale_size_manual(values=c(0.5,0.5,1)) + theme_classic() +xlab(paste(Condition1)) +ylab(paste(Condition2))
  
  Comparison_Data=list("Comparison_table"=Data_test, "scatter_plot"=scatter_plot)
  Comparison_Data
}