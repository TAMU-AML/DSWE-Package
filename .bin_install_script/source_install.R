#################################################################################
#AUTHOR: ABHINAV PRAKASH
#DESCRIPTION: This script can be used to install DSWE package from its source file
##################################################################################

.PACKAGE_VERSION = "1.3.2" #Update this line when version changes
.DEPENDENCIES = list(
  list(name = "Rcpp",version = "1.0.4.6"),
  list(name = "matrixStats",version = "0.55.0"),
  list(name = "FNN",version = "1.1.3"),
  list(name = "KernSmooth",version = "2.23-16"),
  list(name = "mixtools",version = "1.1.0")
)


.FILENAME = paste0("https://github.com/TAMU-AML/DSWE-Package/tarball/v",.PACKAGE_VERSION)  


.packageList = utils::installed.packages()[,1]

for (.i in c(1:length(.DEPENDENCIES))){
  if (.DEPENDENCIES[[.i]]$name %in% .packageList){
    .available_version = utils::packageVersion(.DEPENDENCIES[[.i]]$name)
    if (.available_version < .DEPENDENCIES[[.i]]$version){
      cat(.DEPENDENCIES[[.i]]$name,"version",as.character(.available_version),"is installed, but",.DEPENDENCIES[[.i]]$version,"required",'\n')
      cat("Updating package",.DEPENDENCIES[[.i]]$name,'\n')
      cat("Please select a version higher or equal to",.DEPENDENCIES[[.i]]$version,", if asked",'\n')
      cat("Please select 'compile from source' if asked",'\n')
      install.packages(.DEPENDENCIES[[.i]]$name)
    }
  } else {
    cat(.DEPENDENCIES[[.i]]$name,"is not installed",'\n')
    cat("Installing package",.DEPENDENCIES[[.i]]$name,'\n')
    cat("Please select a version higher or equal to",.DEPENDENCIES[[.i]]$version,", if asked",'\n')
    cat("Please select 'compile from source', if asked",'\n')
    utils::install.packages(.DEPENDENCIES[[.i]]$name)
  }
}

utils::install.packages(.FILENAME,repos=NULL,type = "source")
rm(list = c(".PACKAGE_VERSION",".DEPENDENCIES",".FILENAME",".packageList",".i"))