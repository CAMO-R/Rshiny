rm(list=ls())
#Download R Shiny from github and extract it to a folder "Rshiny". 
#Set the working directory to be the path to the folder "Rshiny".
setwd("path_to_Rshiny/")
shiny::runApp('Rshiny', port=9987, launch.browser=T)


