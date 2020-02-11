##### CENTER AVERAGE AND CORRECTED GENE FITNESS VALUES AROUND THE MODE OF EACH DISTRIBUTION #####

#For each condition, centers the replicate-averaged and mean-corrected fitness values by the mode of the distributio
#This is based on the assumption that most insertions lead to neutral phenotypes, thus the most represented fitness values should be representative of a neutral fitness
#For ease of interpretation, we set up that neutral fitness to 0 by subtracting the mode of a distribution to associated fitness values


Mode_Centered=function(Mean_corrected_averageF,org_locId){
  if(org_locId=="Num"){
    Cond_list=c(unique(as.character(Mean_corrected_averageF$Cdt)))
    Mat_Centered=matrix(0,ncol=ncol(Mean_corrected_averageF))
    colnames(Mat_Centered)=c("locusID","ModeCentered","WeighrtedVar","Cdt")
  
    for (i in 1 : length(Cond_list)){
      Cond=Cond_list[i]
      Data=subset(Mean_corrected_averageF,Mean_corrected_averageF$Cdt==Cond)
      d = density(Data$CorrectedFit,na.rm=TRUE)
      mode = d$x[which.max(d$y)]
    
      Mat_Temp=matrix(nrow=dim(Data)[1],ncol=(ncol(Mean_corrected_averageF)))
      colnames(Mat_Temp)=c("locusID","ModeCentered","WeighrtedVar","Cdt")
      Mat_Temp[,1]=Data[,1]
      Mat_Temp[,2]=Data$CorrectedFit-mode
      Mat_Temp[,3]=Data[,3]
      Mat_Temp[,4]=as.character(Data[,4])
    
      Mat_Centered=rbind(Mat_Centered,Mat_Temp)
    }
  
    Table_Fit_MeanCentered=data.frame("locusId"=as.numeric(Mat_Centered[,1]),"CenteredFit"=as.numeric(Mat_Centered[,2]),
                                    "WeightedVar"=as.numeric(Mat_Centered[,3]),"Cdt"=as.character(Mat_Centered[,4]))
  
    Table_Fit_MeanCentered=Table_Fit_MeanCentered[-1,]
  }
  if(org_locId=="Char"){
    Cond_list=c(unique(as.character(Mean_corrected_averageF$Cdt)))
    Mat_Centered=matrix(0,ncol=ncol(Mean_corrected_averageF))
    colnames(Mat_Centered)=c("locusID","ModeCentered","WeighrtedVar","Cdt")
    
    for (i in 1 : length(Cond_list)){
      Cond=Cond_list[i]
      Data=subset(Mean_corrected_averageF,Mean_corrected_averageF$Cdt==Cond)
      d = density(Data$CorrectedFit,na.rm=TRUE)
      mode = d$x[which.max(d$y)]
      
      Mat_Temp=matrix(nrow=dim(Data)[1],ncol=(ncol(Mean_corrected_averageF)))
      colnames(Mat_Temp)=c("locusID","ModeCentered","WeighrtedVar","Cdt")
      Mat_Temp[,1]=as.character(Data[,1])
      Mat_Temp[,2]=Data$CorrectedFit-mode
      Mat_Temp[,3]=Data[,3]
      Mat_Temp[,4]=as.character(Data[,4])
      
      Mat_Centered=rbind(Mat_Centered,Mat_Temp)
    }
    
    Table_Fit_MeanCentered=data.frame("locusId"=as.character(Mat_Centered[,1]),"CenteredFit"=as.numeric(Mat_Centered[,2]),
                                      "WeightedVar"=as.numeric(Mat_Centered[,3]),"Cdt"=as.character(Mat_Centered[,4]))
    
    Table_Fit_MeanCentered=Table_Fit_MeanCentered[-1,]
  }
  
  Table_Fit_MeanCentered
}