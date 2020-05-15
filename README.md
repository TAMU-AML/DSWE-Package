<center> <h1>DSWE (Data Science for Wind Energy)</h1> </center>


- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#Usage)   
    1. [ComparePCurve](#ComparePCurve)  
    2. [funGP](#funGP)
    3. [KnnPCFit](#knnPCFit)  
    4. [KnnPredict](#KnnPredict)  
    5. [KnnUpdate](#KnnUpdate)  
    6. [AMK](#AMK)  
    7. [CovMatch](#CovMatch)  
- [Details](#details)

# Introduction
This is an R-package implementing some of the data science methods for wind energy applications (DSWE). The current functionalities include creating a multi-dimensional power curve model, performing power curve function comparison, and covariate matching:

Power curve comparison:

* ComparePCurve
* funGP

Predictive modelling functions:

* KnnFit
* KnnPredict
* KnnUpdate
* AMK

Covariate matching function :

* CovMatch

# Installation
The package building relies on certain tool chains in Windows and Mac respectively, as the compiler for C++ code, along with package `remotes`

**Step 1 (Download necessary tool chain):**

Tool chain : [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows, [Apple Command Line Tools](https://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/) and [gFortran](https://mac.r-project.org/tools/) for Mac OS

**Step 2 (Install package remotes in R):**

Install package remotes:
```R
install.packages("remotes")
```

**Step 3 (Build package using remotes):**

```R
remotes::install_github("TAMU-AML/DSWE-Package")
```

You can also specify the version (for example, version 1.3.1) as follows:
```R
remotes::install_github("TAMU-AML/DSWE-Package@v1.3.1")
```

# Usage

The package can be accessed either by attaching the package or using the package name.

*Attaching package and accessing functions*

```R
library(DSWE)

ComparePCurve()
funGP()
KnnFit()
KnnPredict()
KnnUpdate()
AMK()
CovMatch()
```

*Accesing functions without attaching package*
```R
DSWE::ComparePCurve()
DSWE::funGP()
DSWE::KnnFit()
DSWE::KnnPredict()
DSWE::KnnUpdate()
DSWE::AMK()
DSWE::CovMatch()
```

### 1. ComparePCurve
The function can be used to quantify the difference using CovMatch and funGP functions internally.

*Function :*

*ComparePCurve(data, xCol, xCol.circ = NULL, yCol, testCol, testSet = NULL, thrs = 0.2, conflevel = 0.95, gridSize = c(50, 50), powerbins = 15, baseline = 1, limitMemory = T)*

```R
# Preparing the arguments
data1 = read.csv('data1.csv')
data2 = read.csv('data2.csv')
data = list(data1, data2)
xCol = c(1, 3, 6)
xCol.circ = NULL
yCol = 4
testCol = c(1, 3)
testSet = NULL
thrs = 0.2
confLevel = 0.95
gridSize = c(50, 50)
powerbins = 15
baseline = 1
limitMemory = TRUE

# Executing the function
function_comparison = ComparePCurve(data, xCol, xCol.circ, yCol, testCol, testSet, thrs, confLevel, gridSize, powerbins, baseline, limitMemory)
```

### 2. funGP
 The function can be used to perform function comparison using Gaussian process and hypothesis testing

 *Function :*

 *funGP (datalist, xCol, yCol, confLevel = 0.95, testset, limitMemory = TRUE)*

 ```R
 # Preparing the arguments
 datalist = matched_data$matchedData
 xCol = c(1, 3)
 yCol = 4
 confLevel = 0.95
 testset = read.csv('testset.csv')
 limitMemory = TRUE

 # Executing the function
 function_diff = funGP(datalist, xCol, yCol, confLevel, testset, limitMemory)
 ```

### 3. KnnPCFit
The function can be used to model the data using user supplied arguments, a knn model is returned as an end result. It can also be used to get the best feature subset, if subsetSelection is set TRUE

*Function :*

*KnnPCFit(data, xCol, yCol, subsetSelection = FALSE)*

```R
# Preparing the arguments
data = read.csv('data1.csv')
xCol = c(1, 3)
yCol = 4
subsetSelection = FALSE

# Executing the function
knn_model = KnnPCFit(data, xCol, yCol, subsetSelection)
```
### 3. KnnPredict
The function can be used to evaluate a prediction on a new test point using model generated using KnnPCFit

*Function :*

*KnnPredict(knnMdl, testData)*

```R
# Preparing the arguments
knnMdl = knn_model
testData = data[1:100, ]

# Executing the function
prediction = KnnPredict(knnMdl, testData)
```

### 5. KnnUpdate
The function can be used to update the knn model whenever new data in available

*Function :*

*KnnUpdate = function(knnMdl, newData)*

```R
# Preparing the arguments
knnMdl = knn_model
newData = data[500:1000, ]

# Executing the function
knn_newmodel = KnnUpdate(knnMdl, newData)
```

### 6. AMK
The function can be used to model the data by using user supplied arguments. It uses a kernel to assign weights to every training data points, the bandwidth of kernel (bw) can be provided as vector of values or character 'dpi' and 'dpi-gap'. If provided character input, the bandwidths are computed internally

*Function :*

*AMK(trainX, trainY, testX, bw = 'dpi_gap', nMultiCov = 3, fixedCov = c(1, 2), cirCov = NA )*

```R
# Preparing the arguments
data = read.csv('data1.csv')
trainX = data[, c(1, 3)]
trainY = data[, 4]
testX = data[100:200, c(1, 3)]
bw = 'dpi_gap'
nMultiCov = 2
fixedCov = NA
cirCov = NA
# Executing the function
AMK_prediction = AMK(trainX, trainY, testX, bw, nMultiCov, fixedCov, cirCov)
```
**Note :-** In case an error such as a non-finite bandwidth is generated upon adding a covariate, please remove such covariates


### 7. CovMatch
The function can be used to match different data sets. It can only be used to match two different data set at one time. If priority argument is set to FALSE, which is default, the feature columns provided are used in the same order in matching, else computes the covariates matching sequence

*Function :*

*CovMatch(data, xCol, xCol.circ = NULL, thrs = 0.2, priority = FALSE)*

```R
# Preparing the arguments
data1 = read.csv('data1.csv')
data2 = read.csv('data2.csv')

data = list(data1, data2)
xCol = c(1, 3, 6)
xCol.circ = NULL
thrs = c(0.1, 0.1, 0.05)
priority = FALSE

# Executing the function
matched_data = CovMatch(data, xCol, xCol.circ, thrs, priority)
```

**Note :-** Arguments usage detail for each of the functions can be accessed through R documentation using:

```R
help(functionname)
?functionname
```
# Details
For more information on DSWE Package, please access the package documentations. Please feel free to contact the authors.

* Name : Abhinav Prakash, Nitesh Kumar

* Email : abhinavp@tamu.edu, nitesh.kumar@tamu.edu
