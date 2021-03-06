##### COUNTS CORRECTION FILTERING LOW COUNTS< BAD INSERTIONS AND USING A REFERENCE GENE #####

Ref_counts_CorrandNorm=function(Data, ref, cdnb){
  
 # First we want to add a pseudo count of 0.1 to all counts
 for (i in 8:dim(Data)[2]){
   Data[,i]=(Data[,i]+0.1)
 }
 
 # We want to remove any barcode associated with T0 counts <3.1 and f<0.1 and f>0.9
 Data=subset(Data, Data$T0>3.1 & Data$f>0.1 & Data$f<0.9)

 # We calculate the average number of counts per insertion mutant for each individual conditions using a set of provided reference genes 
 ref_data=subset(Data, Data$locusId%in%ref)
 ref_data_counts=ref_data[,8:ncol(ref_data)]
 Scounts_ref=c()

 for (i in 1:dim(ref_data_counts)[2]){
   Scounts_ref[i]=mean(ref_data_counts[,i])
 }

 # We correct the counts for each insertion mutant using the count of the reference gene
 data_counts=Data[,8:ncol(Data)]

 for(i in 1:dim(data_counts)[2]){
   Cor=Scounts_ref[i]
   data_counts[,i]=data_counts[,i]/Cor
  
 }

 Data_Counts_corrected=cbind(Data[1:7],data_counts)
 Data_Counts_corrected

}