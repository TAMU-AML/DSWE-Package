<center> <h1>DSWE (Data Science for Wind Energy)</h1> </center>

# Tutorial and examples

## A tutorial for the package is available on zenodo at this [link](https://zenodo.org/record/6823803#)

# Installation

## For versions 1.5.1 and above (through [CRAN](https://cran.r-project.org/package=DSWE)):
```R
install.packages("DSWE")
```

<span style="color:red"> Note: </span> The package contains C++ code, hence installation using source requires C++ compilers: [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows, [Apple Command Line Tools](https://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/) and [GFortran](https://mac.r-project.org/tools/) for Mac OS. The pre-compiled binaries available through [CRAN](https://cran.r-project.org), the offical package repository for `R`, do not require C++ compilers.


## For versions prior to 1.5.1:
Install using the remotes package by specify the version (for example, version 1.3.1) as follows:
```R
remotes::install_github("TAMU-AML/DSWE-Package@v1.3.1")
```
