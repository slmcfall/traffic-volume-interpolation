


# functions

# sets working directory
# made so I don't have to copy paste a setwd() funciton all the time
set_directory <- function(computer_name){
  if (computer_name == "windowsCGA"){
    setwd("W:/slmcfall/traffic_from_desktop/ID/Pri_Sec_combined2/")
    
  } else if (computer_name == "linuxCGA"){
    setwd("/media/sean/Data/Research/trafficVolume/ID/Pri_Sec_combined2/")
    
  } else if (computer_name == "SWEMCGA18"){
    setwd("W:/slmcfall/traffic_from_desktop/ID/Pri_Sec_combined2/")
  }
}

