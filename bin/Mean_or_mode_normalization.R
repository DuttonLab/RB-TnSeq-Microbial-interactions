##### NORMALIZATION 2 BY CENTERING AROUND THE MEAN OR THE MODE #####

# This script contains two functions: one to normalize location-normalized fitness values around the mean of the distribution and one to centered the location-normalized fitness values around the mode of the distribution

# Last step of normalization is, per condition, to center each fitness value around either the mode of the distribution or the mean
#Function to normalize around the mean
Norm_around_mean=function(Gene_Fit_chromNorm,cdnbF){

  Data_centered_mean=Gene_Fit_chromNorm
  mean_vec=vector(mode="numeric", length=cdnbF)
  names(mean_vec)=colnames(Data_centered_mean[,5:ncol(Data_centered_mean)])

  for (i in 5:ncol(Data_centered_mean)){
    Data_centered_mean[,i]=as.numeric(Data_centered_mean[,i])
    meanF=mean(Data_centered_mean[,i])
    Data_centered_mean[,i]=(Data_centered_mean[,i]-meanF)
    mean_vec[i-4]=meanF
  }
  Data_output=list(Data_centered_mean,mean_vec)
  Data_output
}

# Function to center around the mode
Norm_around_mode=function(Gene_Fit_chromNorm){

  Data_centered_mode=Gene_Fit_chromNorm

  for (i in 5:ncol(Data_centered_mode)){
    Data_centered_mode[,i]=as.numeric(Data_centered_mode[,i])
    d = density(Data_centered_mode[,i])
    mode = d$x[which.max(d$y)]
    Data_centered_mode[,i]=(Data_centered_mode[,i]-mode)
  }
  Data_centered_mode
}
