##### DATA PREPARATION FOR VISUALIZATION OF COUNTS PER 10KB #####

Data_prep_viz10KB=function(Data,cdnb, scaffoldX){

  Data=subset(Data, Data$scaffold==scaffoldX)
# How many 10kB section is there in the library
  position_list=c(Data$pos)
  frag=round((max(position_list)-min(position_list))/10000) #number of fragment of 10000kb ==> used for the for loop

#How many conditions and names     # we determine the number of conditions (including the T0) in your study
  cd_names=c(colnames(Data[,8:ncol(Data)]))   # create a vector with the names of the conditions

#We create an empty matrix to further store number of reads per 10kB in each condition
  Mat_10KbC=matrix(0, nrow=frag, ncol=cdnb)
  colnames(Mat_10KbC)=cd_names

  BC_per10kb=c()
#We fill the matrix. We counts the number of reads per 10kB per condition by chunking the Data table every 10Kb based on the pos column

  for (i in 1:frag){
    Dat_kb=subset(Data, Data$pos<=10000)
    Dat_kb=Dat_kb[,-1:-7]
    len=dim(Dat_kb)[1]
    for (j in 1 : cdnb){
      Mat_10KbC[i,j]=sum(Dat_kb[,j])
    
    }
    BC_per10kb[i]=len
    Data=Data[-1:-len,]
    Data$pos=(Data$pos-10000)
  }

  Count_10kb_data=as.data.frame(Mat_10KbC)
  Count_10kb_data$Rank=c(1:frag)
  
  Count_10kb_data

}

