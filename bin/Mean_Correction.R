##### CORRECTION OF THE AVERAGE GENE FITNESS BY THE MEAN OF THE REPLICATE MEANS #####

#This function allows, for each condition, to correct the replicate-averaged fitness values by the Mean of the replicates mean fitness of the distribution (saved in the Rdata file) 


Mean_Correction=function(Average_fitness, mean1,mean2,mean3,org_locId){
  if(org_locId=="Num"){
    Cond_list=c(unique(as.character(Average_fitness$Cdt)))
    Mat_Corrected=matrix(0,ncol=(ncol(Average_fitness)+1))
    colnames(Mat_Corrected)=c(colnames(Average_fitness),"WeightedFitCor")
  
    for (i in 1 : length(Cond_list)){
      Cond=Cond_list[i]
      Data=subset(Average_fitness, Average_fitness$Cdt==Cond)
      m1=mean1[Cond]
      m2=mean2[Cond]
      m3=mean3[Cond]
      MeanR=mean(c(m1,m2,m3))
    
      Mat_Temp=matrix(nrow=dim(Data)[1],ncol=(ncol(Average_fitness)+1))
      colnames(Mat_Temp)=c(colnames(Average_fitness),"WeightedFitCor")
      Mat_Temp[,1]=Data[,1]
      Mat_Temp[,2]=Data[,2]
      Mat_Temp[,3]=Data[,3]
      Mat_Temp[,4]=as.character(Data[,4])
      Mat_Temp[,5]=Data$WeightedFit+MeanR
    
      Mat_Corrected=rbind(Mat_Corrected,Mat_Temp)
    }
  
    Table_Fit_MeanCorrected=data.frame("locusId"=as.numeric(Mat_Corrected[,1]),"CorrectedFit"=as.numeric(Mat_Corrected[,5]),
                                     "WeightedVar"=as.numeric(Mat_Corrected[,3]),"Cdt"=as.character((Mat_Corrected[,4])))
  
    Table_Fit_MeanCorrected=Table_Fit_MeanCorrected[-1,]
  }
  
  if(org_locId=="Char"){
    Cond_list=c(unique(as.character(Average_fitness$Cdt)))
    Mat_Corrected=matrix(0,ncol=(ncol(Average_fitness)+1))
    colnames(Mat_Corrected)=c(colnames(Average_fitness),"WeightedFitCor")
    
    for (i in 1 : length(Cond_list)){
      Cond=Cond_list[i]
      Data=subset(Average_fitness, Average_fitness$Cdt==Cond)
      m1=mean1[Cond]
      m2=mean2[Cond]
      m3=mean3[Cond]
      MeanR=mean(c(m1,m2,m3))
      
      Mat_Temp=matrix(nrow=dim(Data)[1],ncol=(ncol(Average_fitness)+1))
      colnames(Mat_Temp)=c(colnames(Average_fitness),"WeightedFitCor")
      Mat_Temp[,1]=as.character(Data[,1])
      Mat_Temp[,2]=Data[,2]
      Mat_Temp[,3]=Data[,3]
      Mat_Temp[,4]=as.character(Data[,4])
      Mat_Temp[,5]=Data$WeightedFit+MeanR
      
      Mat_Corrected=rbind(Mat_Corrected,Mat_Temp)
    }
    
    Table_Fit_MeanCorrected=data.frame("locusId"=as.character(Mat_Corrected[,1]),"CorrectedFit"=as.numeric(Mat_Corrected[,5]),
                                       "WeightedVar"=as.numeric(Mat_Corrected[,3]),"Cdt"=as.character((Mat_Corrected[,4])))
    
    Table_Fit_MeanCorrected=Table_Fit_MeanCorrected[-1,]
  }

  Table_Fit_MeanCorrected
}