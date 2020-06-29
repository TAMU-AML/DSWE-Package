#################################################################################
#AUTHOR: ABHINAV PRAKASH
#DESCRIPTION: This script can be used to install DSWE package from its source file
##################################################################################

.PACKAGE_VERSION = "1.3.4" #Update this line when version changes
.requiredRversion = "3.5.0"
.DEPENDENCIES = list(
  list(name = "Rcpp",version = "1.0.4.6"),
  list(name = "matrixStats",version = "0.55.0"),
  list(name = "FNN",version = "1.1.3"),
  list(name = "KernSmooth",version = "2.23-16"),
  list(name = "mixtools",version = "1.1.0"),
  list(name = "RcppArmadillo",version = "0.9.870.2.0")
)


.Rversion = paste0(R.Version()$major,".",R.Version()$minor)
if (.Rversion < .requiredRversion){
  message("Aborting Installation ... ")
  stop("R version 3.5.0 or higher is required for DSWE package", call. = FALSE)
}
message("Proceeding Installation ... ")
.packageList = utils::installed.packages()[,1]
for (.i in c(1:length(.DEPENDENCIES))){
  if (.DEPENDENCIES[[.i]]$name %in% .packageList){
    .available_version = utils::packageVersion(.DEPENDENCIES[[.i]]$name)
    if (.available_version < .DEPENDENCIES[[.i]]$version){
      message(.DEPENDENCIES[[.i]]$name,"version",as.character(.available_version),"is installed, but",.DEPENDENCIES[[.i]]$version,"required")
      message("Updating package",.DEPENDENCIES[[.i]]$name)
      message("Please select a version higher or equal to",.DEPENDENCIES[[.i]]$version,", if asked")
      message("Please select 'compile from source' if asked")
      install.packages(.DEPENDENCIES[[.i]]$name)
    }
  } else {
    message(.DEPENDENCIES[[.i]]$name,"is not installed")
    message("Installing package",.DEPENDENCIES[[.i]]$name)
    message("Please select a version higher or equal to",.DEPENDENCIES[[.i]]$version,", if asked")
    message("Please select 'compile from source', if asked")
    utils::install.packages(.DEPENDENCIES[[.i]]$name)
  }
}
.downloadURL = paste0("https://github.com/TAMU-AML/DSWE-Package/archive/v",.PACKAGE_VERSION,".tar.gz") 
.rootdir = tempfile(pattern = "DSWE")
dir.create(.rootdir)
.destfile = tempfile(pattern = "DSWE", fileext = ".tar.gz", tmpdir = .rootdir)
download.file(url = .downloadURL, destfile = .destfile)
untar(.destfile, exdir = .rootdir )
.folder = list.dirs(.rootdir, full.names = T, recursive = F)
utils::install.packages(.folder,repos=NULL,type = "source")
rm(list = c(".PACKAGE_VERSION",".DEPENDENCIES",".downloadURL",".rootdir",".destfile",".folder",".packageList",".i"))

