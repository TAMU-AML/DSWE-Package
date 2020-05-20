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
  list(name = "mixtools",version = "1.1.0"),
  list(name = "RcppArmadillo",version = "0.9.870.2.0")
)

message("You are about to install DSWE package. Please make sure that you are using R version 3.5.0 or higher")
message("Press enter to continue, or any other key and enter to abort:")
.action = readline()
if (.action != ""){
  rm(".PACKAGE_VERSION",".DEPENDENCIES",".action")
  message("Aborting Installation ... ")
  stop("Aborted by user", call. = FALSE)
} else {
  message("Proceeding Installation ... ")
  if (.Platform$OS.type == "windows" ){
      .Platform$file.sep = "\\\\"
    .foundRtools = grepl("Rtools",Sys.getenv("PATH"))
    if (!.foundRtools) {
      message("We could not find Rtools in 'PATH', Please select an appropriate action:")
      cat("1:","Abort installation. I don't have Rtools installed.",'\n')
      cat("2:","Add Rtools to PATH",'\n')
      cat("3:","Continue anyway. I am confident that I have Rtools in PATH",'\n')
      .action = readline("Enter a number indicating your response: ")
      if (!(as.integer(.action) %in% c(1,2,3))){
        stop("Invalid option, aborting installation ... ", call. = FALSE)
      } else if (as.integer(.action) == 1){
         rm(".PACKAGE_VERSION",".DEPENDENCIES",".action",".foundRtools")
         message("Aborting Installation ... ")
         stop("Aborted by user", call. = FALSE)
      } else if (as.integer(.action) == 2 || as.integer(.action) == 3){
        if (as.integer(.action) == 2){
          .pathRtools = readline("Enter the path to Rtools, C:\\Rtools\\bin, for example: ")
          Sys.setenv("PATH" = paste(.pathRtools,Sys.getenv("PATH"),sep = ";"))
          cat("Rtools added to PATH",'\n')
        }
      } 
    }
  }
  .packageList = utils::installed.packages()[,1]
  for (.i in c(1:length(.DEPENDENCIES))){
    if (.DEPENDENCIES[[.i]]$name %in% .packageList){
      .available_version = utils::packageVersion(.DEPENDENCIES[[.i]]$name)
      if (.available_version < .DEPENDENCIES[[.i]]$version){
        message(.DEPENDENCIES[[.i]]$name,"version",as.character(.available_version),"is installed, but",.DEPENDENCIES[[.i]]$version,"required",'\n')
        message("Updating package",.DEPENDENCIES[[.i]]$name,'\n')
        message("Please select a version higher or equal to",.DEPENDENCIES[[.i]]$version,", if asked",'\n')
        message("Please select 'compile from source' if asked",'\n')
        install.packages(.DEPENDENCIES[[.i]]$name)
      }
    } else {
      message(.DEPENDENCIES[[.i]]$name,"is not installed",'\n')
      message("Installing package",.DEPENDENCIES[[.i]]$name,'\n')
      message("Please select a version higher or equal to",.DEPENDENCIES[[.i]]$version,", if asked",'\n')
      message("Please select 'compile from source', if asked",'\n')
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
  rm(list = c(".PACKAGE_VERSION",".DEPENDENCIES",".downloadURL",".rootdir",".destfile",".folder",".action",".packageList",".i"))
}
