#################################################################################
#AUTHOR: ABHINAV PRAKASH
#DESCRIPTION: This script can be used to install DSWE package from its binary file
##################################################################################

.PACKAGE_VERSION = "1.3.4" #Update this line when version changes
.DEPENDENCIES = list(
  list(name = "Rcpp",version = "1.0.4.6"),
  list(name = "matrixStats",version = "0.55.0"),
  list(name = "FNN",version = "1.1.3"),
  list(name = "KernSmooth",version = "2.23-16"),
  list(name = "mixtools",version = "1.1.0")
)

if (Sys.info()["sysname"]=="Windows"){
  .FILENAME = paste0("https://github.com/TAMU-AML/DSWE-Package/releases/download/v",.PACKAGE_VERSION,"/DSWE_",.PACKAGE_VERSION,".zip")  
} else {
  .FILENAME = paste0("https://github.com/TAMU-AML/DSWE-Package/releases/download/v",.PACKAGE_VERSION,"/DSWE_",.PACKAGE_VERSION,".tgz")  
}

.packageList = utils::installed.packages()[,1]

for (.i in c(1:length(.DEPENDENCIES))){
  if (.DEPENDENCIES[[.i]]$name %in% .packageList){
    .available_version = utils::packageVersion(.DEPENDENCIES[[.i]]$name)
    if (.available_version < .DEPENDENCIES[[.i]]$version){
      cat(.DEPENDENCIES[[.i]]$name,"version",as.character(.available_version),"is installed, but",.DEPENDENCIES[[.i]]$version,"required",'\n')
      cat("Updating package",.DEPENDENCIES[[.i]]$name,'\n')
      install.packages(.DEPENDENCIES[[.i]]$name)
    }
  } else {
    cat(.DEPENDENCIES[[.i]]$name,"is not installed",'\n')
    cat("Installing package",.DEPENDENCIES[[.i]]$name,'\n')
    utils::install.packages(.DEPENDENCIES[[.i]]$name)
  }
}

utils::install.packages(.FILENAME,repos=NULL)

rm(list = c(".PACKAGE_VERSION",".DEPENDENCIES",".FILENAME",".packageList",".i"))
