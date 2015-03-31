#validationMethods.R

sampling<-vector("list",length(ghcn.subsets))
sampling_station_id<-vector("list",length(ghcn.subsets))
for(i in 1:length(ghcn.subsets)){
  n<-nrow(ghcn.subsets[[i]])
  prop<-(sampling_dat$prop[i])/100
  ns<-n-round(n*prop)   #Create a sample from the data frame with 70% of the rows
  nv<-n-ns              #create a sample for validation with prop of the rows
  ind.training <- sample(nrow(ghcn.subsets[[i]]), size=ns, replace=FALSE) #This selects the index position for 70% of the rows taken randomly
  ind.testing <- setdiff(1:nrow(ghcn.subsets[[i]]), ind.training)
  #Find the corresponding 
  data_sampled<-ghcn.subsets[[i]][ind.training,] #selected the randomly sampled stations
  station_id.training<-data_sampled$station     #selected id for the randomly sampled stations (115)
  #Save the information
  sampling[[i]]<-ind.training #index of training sample from data.frame
  sampling_station_id[[i]]<- station_id.training #station ID for traning samples
}