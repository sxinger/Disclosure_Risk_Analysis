##---------------------------helper functions--------------------------------------##
## install (if needed) and require packages
require_libraries<-function(package_list){
  #install missing packages
  install_pkg<-as.data.frame(installed.packages())
  new_packages<-package_list[!(package_list %in% install_pkg[which(install_pkg$LibPath==.libPaths()[1]),"Package"])]
  if(length(new_packages)>0){
    install.packages(new_packages,lib=.libPaths()[1],repos = "http://cran.us.r-project.org")
  }
  
  for (lib in package_list) {
    library(lib, character.only=TRUE,lib.loc=.libPaths()[1])
    cat("\n", lib, " loaded.", sep="")
  }
}

## procecutor attacker model risk functions
#https://github.com/sdcTools/sdcMicro/blob/9d4b05193ec5c4db0745bacb6ac470d6b358a363/R/modRisk.R#L138
#1. estimates the number of sample uniques that are population unique
risk1 <- function(l, p) {
  v=(1 - p) * l
  exp(-v)
} 

#2. estimates the number of correct matches of sample uniques
risk2 <- function(l, p) {
  v=(1 - p) * l
  (1 - exp(-v))/v
}