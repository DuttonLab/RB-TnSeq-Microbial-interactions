##### NORMALIZATION BY SMOOTHED MEDIAN #####


# If your species have chromosome + plasmid, the correction must be done for each scaffold independantly
# To get the smoothed median we use a window of 251 genes, if one scaffold is less than 251 genes but more than 10 genes, we correct by the median as described in Wetmore et al., 2015
# If the scaffold is less than 10 genes, we do not correct

Loc_smoothmed_norm=function(Table_GeneFitness_MEAN, genes.tab, cdnbF,org_locId){

# First we detrmine the number of scaffolds
  
  Data_norm1=Table_GeneFitness_MEAN
  
  genes.tab=genes.tab[order(genes.tab$locusId),]
  genes.tab=subset(genes.tab, genes.tab$locusId%in%Data_norm1$locusId)
  
  Data_norm1$scaffoldId=genes.tab$scaffoldId

  scaffoldId=unique(Data_norm1$scaffoldId)
  nb_scaffold=length(scaffoldId)

#We set up an empty matrix we will fill with the normalized values

  Mat_norm=matrix(0,nrow=1, ncol=dim(Data_norm1)[2])
  cd_names=colnames(Data_norm1[2:(ncol(Data_norm1)-3)])
  colnames(Mat_norm)=c("locusId","sysName","begin","scaffoldId",cd_names)

  window=251  # Have to use a window with an odd number for smoothed median
  c=cdnbF  #number of condition (T0 removed)

  for (i in 1:nb_scaffold){
    scaf=scaffoldId[[i]]
    Dat=subset(Data_norm1,Data_norm1$scaffoldId==scaf)
    gene_sca=dim(Dat)[1]
    Mat_temp=matrix(0,nrow=dim(Dat)[1], ncol=dim(Dat)[2])
    colnames(Mat_temp)=c("locusId","sysName","begin","scaffoldId",cd_names)
  
    Mat_temp[,1]=Dat$locusId
    if (org_locId=="Char"){
      Mat_temp[,1]=as.character(Dat$locusId)
    }
    Mat_temp[,2]=Dat$sysName
    Mat_temp[,3]=Dat$begin
    Mat_temp[,4]=Dat$scaffoldId
   
    if(gene_sca>=window){
      for(j in 2:(2+c-1)){
        values=as.numeric(Dat[,j])
        names(values)=Dat[,1]
    
        med=runmed(values,window,endrule="constant")
    
        values_corr=values-med
    
        Mat_temp[,(3+j)]=values_corr
      }
    }
    
    if (gene_sca>=10 & gene_sca<251){
      for(j in 2:(2+c-1)){
        values=as.numeric(Dat[,j])
        names(values)=Dat[,1]
        
        med=median(values)
        
        values_corr=values-med
        
        Mat_temp[,(3+j)]=values_corr
      } 
    }
    if (gene_sca<10){
      for(j in 2:(2+c-1)){
        values=as.numeric(Dat[,j])
        names(values)=Dat[,1]
        Mat_temp[,(3+j)]=values
      } 
    }
    Mat_norm=rbind(Mat_norm,Mat_temp)
  }

  Gene_Fit_chromNorm=as.data.frame(Mat_norm[-1,],stringsAsFactors=FALSE)
  Gene_Fit_chromNorm

}
