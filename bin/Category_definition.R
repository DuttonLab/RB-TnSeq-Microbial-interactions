### Determines the category of a gene in a 2condition comparison ###

## Possible categories: "Not_testd": the test for equal variance was not met, so no comparison has been performed for that gene
                       #"Not_sig": no significant difference has been observed for that gene between the 2 tested conditions
                      #  "Sig": that gene has a significant fitness difference between the conditions

Category_definition=function(Data_test,alpha, multi){
  
  if (multi==0){
    Cutoff=abs(qt(alpha/2,4))
  
    Data_test$Category="Not_sig"
  
    for (i in 1:dim(Data_test)[1]){
      if (is.na(Data_test[i,17])==TRUE){
        Data_test[i,20]="Not_tested"
      }
      if (is.na(Data_test[i,17])==FALSE){
       Data_test[i,20]="Not_Sig"
      
        if (abs(Data_test[i,17])>=Cutoff){
          Data_test[i,20]="Sig"
        }
      }
    }
  }
  
  if (multi==1){
    
    Data_test$Category="Not_sig"
    Cutoff=alpha
    
    for (i in 1:dim(Data_test)[1]){
      if (is.na(Data_test[i,19])==TRUE){
        Data_test[i,20]="Not_tested"
      }
      if (is.na(Data_test[i,19])==FALSE){
        Data_test[i,20]="Not_Sig"
        
        if (abs(Data_test[i,19])<=Cutoff){
          Data_test[i,20]="Sig"
        }
      }
    }
  }
  Data_test
}