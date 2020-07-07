##### CALCULATION OF FITNESS VALUES CORRELATION BETWEEEN REPLICATES #####

#This function: (i) calculates for each condition, the gene fitness values correlation between replicates, looking at each possible combination of replicates (Pearson, Spearman and Lin)
#               (ii) produces a pdf file with the correlation plots for each correlation analysis performed
#               (iii) generates a table containing all the correlation coefficients (Pearson (P_Rsquared), Spearman (S_Rsquared) and Lin (L_Rsquared))


Correlation_Rep=function(AllReplicate){
  
  Combi=c(unique(AllReplicate$Rep))
  Combi_mat=combn(Combi,2)
  
  Condition=c(unique(as.character(AllReplicate$Cdt)))
  
  Mat_Storage_Cor=matrix(0,nrow=0,ncol=5)
  colnames(Mat_Storage_Cor)=c("Cond","Comp","Pearson","Spearman","Lins")
  
  pdf("Correlation_plots_fit.pdf")
  for(i in 1:length(Condition)){
    Cond=Condition[i]
    Data_Cond=subset(AllReplicate, AllReplicate$Cdt==Cond)
    Mat_Temp=matrix(0,nrow=dim(Combi_mat)[2],ncol=5)
    colnames(Mat_Temp)=c("Cond","Comp","Pearson","Spearman","Lins")
    for(j in 1:dim(Combi_mat)[2]){
      Comp_set=Combi_mat[,j]
      Set1=subset(Data_Cond, Data_Cond$Rep==Comp_set[1])
      Set2=subset(Data_Cond, Data_Cond$Rep==Comp_set[2])
      CorP_val=cor(Set1$NormFitness, Set2$NormFitness, method=c("pearson"))
      CorS_val=cor(Set1$NormFitness, Set2$NormFitness, method=c("spearman"))
      CorL_val=CCC(Set1$NormFitness, Set2$NormFitness)
      
      Data_plot=left_join(Set1,Set2,by=c("locusId"))
      
      plotcor=ggplot(Data_plot, aes(x=NormFitness.x, y=NormFitness.y )) + geom_point(size=0.5, shape=20)+
        geom_abline(intercept = 0, slope = 1, col="red", linetype='dashed') + theme_light()+
        geom_text(x=0, y=min(Data_plot$NormFitness.y), label=paste("R2 = ", CorP_val), col="red")+
        labs(y = paste(Comp_set[2]), x=paste(Comp_set[1])) + ggtitle(paste(Cond,":",Comp_set[2],"vs",Comp_set[1]))
      
      print(plotcor)
      
      Mat_Temp[,1]=paste(Cond)
      Mat_Temp[j,2]=paste(Comp_set[2],"vs",Comp_set[1])
      Mat_Temp[j,3]=CorP_val
      Mat_Temp[j,4]=CorS_val
      Mat_Temp[j,5]=CorL_val$rho.c$est
      
      
    }
    
    Mat_Storage_Cor=rbind(Mat_Storage_Cor,Mat_Temp)
    Correlation_table=data_frame("Cond"=as.character(Mat_Storage_Cor[,1]),
                                 "Comp"=as.character(Mat_Storage_Cor[,2]),
                                 "P_Rsquared"=as.numeric(Mat_Storage_Cor[,3]),
                                 "S_Rsquared"=as.numeric(Mat_Storage_Cor[,4]),
                                 "L_Rsquared"=as.numeric(Mat_Storage_Cor[,5]))
  }
  dev.off()
  Correlation_table
}