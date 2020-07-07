##### CALCULATION OF THE WEIGHTED AVERAGE BY THE INVERSE OF VARIANCE OF GENE FITNESS  #####

#This function: (i) calculates for each condition and each gene the average gene fitness as the weighted fitness across repliactes by invesre of variance
#               (ii) calculates the associated variance or squared standard error (inverse of the sum of inverse variance)
#               (iii) calculates the associated standard deviation

Weighted_average=function(AllReplicate, org_locId){
  #First "for loop" generates the values used for the weighted average calculation:
  #inV= inverse of the variance
  #fitV= fitness value divided by the variance
  #as some genes don't have a variance value we use a variance of 1 to calculate for these genes the usual arithmetic mean
  AllReplicate_plus=AllReplicate
  AllReplicate_plus$invV=1
  AllReplicate_plus$fitV=AllReplicate_plus$NormFitness
  
  for (i in 1:dim(AllReplicate_plus)[1]){
    NA_test=is.na(AllReplicate_plus[i,7])
    if (NA_test[[1]]!=TRUE & AllReplicate[i,7]!=0){
      AllReplicate_plus[i,9]=(1/AllReplicate_plus[i,7])
      AllReplicate_plus[i,10]=(AllReplicate_plus[i,6]/AllReplicate_plus[i,7])
    }
  }
  
  # Now we calculate the weighted fitness value (weightedFit), the asociated variance (weightedVar) and the associated standard deviation (weightedSdev)
  # Do it per condition and per gene ==> nested "for loop"
  
  Mat_weigthed=matrix(0,ncol=5)
  colnames(Mat_weigthed)=c("locusId","Cond","weightedFit","weightedVar","weightedSdev")
  Cond_list=c(unique(as.character(AllReplicate$Cdt)))
  locus_list=c(unique(AllReplicate$locusId))
  
  for(i in 1:length(Cond_list)){
    Cond=Cond_list[i]
    Data=subset(AllReplicate_plus, AllReplicate_plus$Cdt==Cond)
    Mat_Tem=matrix(0,nrow=length(locus_list),ncol=5)
    colnames(Mat_Tem)=c("locusId","Cond","weightedFit","weightedVar","weightedSdev")
    Mat_Tem[,2]=Cond
    for(j in 1:length(locus_list)){
      Data2=subset(Data, Data$locusId==locus_list[j])
      Wt_Fit=(sum(Data2$fitV))/(sum(Data2$invV)) #calculates the weighted mean
      
      if(sum(Data2$invV)!=3){  #The next 2 if loops calculates the variance as the squared standard error
        Wt_Var=1/(sum(Data2$invV)) 
      }
      
      if (sum(Data2$invV)==3){
        Wt_var=var(Data2$fitV) 
      }
      
      if(sum(Data2$invV)!=3){ # The next 2 if loops calculates the standard deviation (Note: this is the equation for 3 replicates 2/3 ==> (N-1)/N)
        Num=sum((Data2$NormFitness-Wt_Fit)^2*Data2$invV)
        Den=2/3*(sum(Data2$invV))
        Wt_sdev=sqrt(Num/Den)
      }
      
      if (sum(Data2$invV)==3){
        Wt_sdev=sd(Data2$fitV) 
      }
      
      Mat_Tem[j,1]=locus_list[j]
      Mat_Tem[j,3]=Wt_Fit
      Mat_Tem[j,4]=Wt_Var
      Mat_Tem[j,5]=Wt_sdev
    }
    Mat_weigthed=rbind(Mat_weigthed,Mat_Tem)
    
  }
  if(org_locId=="Num"){
    Table_weighted_average=data.frame("locusId"=as.numeric(Mat_weigthed[,1]),"WeightedFit"=as.numeric(Mat_weigthed[,3]),"WeightedVar"=as.numeric(Mat_weigthed[,4]),
                                      "Weightedsdev"=as.numeric(Mat_weigthed[,5]),"Cdt"=as.character(Mat_weigthed[,2]))
    Table_weighted_average=Table_weighted_average[-1,]
  }
  
  if(org_locId=="Char"){
    Table_weighted_average=data.frame("locusId"=as.character(Mat_weigthed[,1]),"WeightedFit"=as.numeric(Mat_weigthed[,3]),"WeightedVar"=as.numeric(Mat_weigthed[,4]),
                                      "Weightedsdev"=as.numeric(Mat_weigthed[,5]),"Cdt"=as.character(Mat_weigthed[,2]))
    Table_weighted_average=Table_weighted_average[-1,]
  }
  
  Table_weighted_average
}