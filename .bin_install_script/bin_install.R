#################################################################################
#AUTHOR: ABHINAV PRAKASH
#DESCRIPTION: This script can be used to install DSWE package from its binary file
##################################################################################

PACKAGE_VERSION = "1.3.1" #Update this line when version changes
DEPENDENCIES = list(
  c(name = "Rcpp",version = "1.0.4.6"),
  c(name = "matrixStats",version = "0.55.0"),
  c(name = "FNN",version = "1.1.3"),
  c(name = "KernSmooth",version = "2.23-16"),
  c(name = "mixtools",version = "1.1.0")
)

if (Sys.info()["sysname"]=="Windows"){
  FILENAME = paste0("https://github.com/TAMU-AML/DSWE-Package/releases/download/v",PACKAGE_VERSION,"/DSWE_",PACKAGE_VERSION,".zip")  
} else {
  FILENAME = paste0("https://github.com/TAMU-AML/DSWE-Package/releases/download/v",PACKAGE_VERSION,"/DSWE_",PACKAGE_VERSION,".tgz")  
}


for (i in 1:length(DEPENDENCIES)){
  package_available = require(DEPENDENCIES[[i]]["name"], quietly = TRUE)
  if(package_available){
    available_version = utils::packageVersion(DEPENDENCIES[[i]]["name"])
    if (available_version < DEPENDENCIES[[i]]["version"]){
      install.packages(DEPENDENCIES[[i]]["name"])
    }
  } else {
    install.packages(DEPENDENCIES[[i]]["name"])
  }
}

install.packages(FILENAME,repos=NULL)
