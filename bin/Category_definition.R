### Determines the category of a gene in a 2condition comparison ###

## Possible categories: "Not_testd": the test for equal variance was not met, so no comparison has been performed for that gene
                       #"Not_sig": no significant difference has been observed for that gene between the 2 tested conditions
                      #  "Sig": that gene has a significant fitness difference between the conditions

Category_definition=function(Data_test,alpha){
  
  Cutoff=abs(qt(alpha/2,4))
  
  Data_test$Category="Not_sig"
  
  for (i in 1:dim(Data_test)[1]){
    if (is.na(Data_test[i,15])==TRUE){
      Data_test[i,16]="Not_tested"
    }
    if (is.na(Data_test[i,15])==FALSE){
      Data_test[i,16]="Not_Sig"
      
      if (abs(Data_test[i,15])>=Cutoff){
        Data_test[i,16]="Sig"
      }
    }
  }
  Data_test
}