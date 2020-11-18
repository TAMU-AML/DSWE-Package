<center> <h1>DSWE (Data Science for Wind Energy)</h1> </center>

 ### <span style="color:red"> Note: </span> A graphical installation and usage guide is available on TAMU Advanced Metrology Lab's website at this [link](http://11ekj91tjuht16uu5226632v-wpengine.netdna-ssl.com/wp-content/uploads/sites/164/2020/05/DSWE_HELP.pdf).</h2>

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#Usage)   
    1. [ComparePCurve](#ComparePCurve) 
    2. [ComputeWeightedDifference](#ComputeWeightedDifference) 
    3. [funGP](#funGP)
    4. [tempGP](#tempGP)
    5. [KnnPCFit](#KnnPCFit)  
    6. [KnnPredict](#KnnPredict)  
    7. [KnnUpdate](#KnnUpdate)  
    8. [AMK](#AMK)  
	9. [BayesTreePCFit](#BayesTreePCFit)
	10. [SplinePCFit](#SplinePCFit)
	11. [SvmPCFit](#SvmPCFit)
    12. [CovMatch](#CovMatch)  
- [Details](#details)

# Introduction
This is an R-package implementing some of the data science methods for wind energy applications (DSWE). The current functionalities include creating a multi-dimensional power curve model, performing power curve function comparison, and covariate matching:

Power curve comparison:

* ComparePCurve
* ComputeWeightedDifference
* funGP


Predictive modelling functions:

* tempGP
* KnnPCFit
* KnnPredict
* KnnUpdate
* AMK
* BayesTreePCFit
* SplinePCFit
* SvmPCFit


Covariate matching function :

* CovMatch

# Installation
The package building relies on certain tool chains in Windows and Mac respectively, as the compiler for C++ code, along with package `remotes`

**Step 1 (Download necessary tool chain):**

Tool chain : [Rtools](https://cran.r-project.org/bin/windows/Rtools/) for Windows, [Apple Command Line Tools](https://osxdaily.com/2014/02/12/install-command-line-tools-mac-os-x/) and [GFortran](https://mac.r-project.org/tools/) for Mac OS

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
ComputeWeightedDifference()
funGP()
tempGP()
KnnPCFit()
KnnPredict()
KnnUpdate()
AMK()
BayesTreePCFit()
SplinePCFit()
CovMatch()
```

*Accesing functions without attaching package*
```R
DSWE::ComparePCurve()
DSWE::ComputeWeightedDifference()
DSWE::funGP()
DSWE::tempGP()
DSWE::KnnPCFit()
DSWE::KnnPredict()
DSWE::KnnUpdate()
DSWE::AMK()
DSWE::BayesTreePCFit
DSWE::SplinePCFit
DSWE::CovMatch()
```

### 1. ComparePCurve
The function can be used to quantify the difference using CovMatch and funGP functions internally.

*Function :*

*ComparePCurve(data, xCol, xCol.circ = NULL, yCol, testCol, testSet = NULL, thrs = 0.2, conflevel = 0.95, gridSize = c(50, 50), powerbins = 15, baseline = 1, limitMemory = T, opt_method = 'L-BFGS-B', sampleSize = list(optimSize = 500, bandSize = 5000), rngSeed = 1)*

```R
# Preparing the arguments
data1 = data1[1:100, ]
data2 = data2[1:100, ]
data = list(data1, data2)
xCol = c(2, 4)
xCol.circ = NULL
yCol = 7
testCol = c(2, 4)
testSet = NULL
thrs = 0.2
confLevel = 0.95
gridSize = c(10, 10)
powerbins = 15
baseline = 1
limitMemory = TRUE
opt_method = 'L-BFGS-B'
sampleSize = list(optimSize = 500, bandSize = 5000) 
rngSeed = 1 

# Executing the function
function_comparison = ComparePCurve(data, xCol, xCol.circ, yCol, testCol, testSet, thrs, confLevel, gridSize, powerbins, baseline, limitMemory, opt_method, sampleSize, rngSeed)
```

### 2. ComputeWeightedDifference
This function is used to compute the weighted difference between power curves using user-specified weights.
*Function :*

*ComputeWeightedDifference(muDiff, weights, base, statDiff = FALSE, confBand = NULL)*


 ```R
 # Construct a desired testset and calculate weights based on any probability distribution
 ws_test = seq(3,15,length.out = 10)  #generate 10 grid points for wind speed.
 density_test = seq(1.12, 1.16, length.out = 10) #generate 10 grid points for temperature.
 
 #Combine ws_test and temp_test to create a 50 by 50 mesh grid.
 testset = expand.grid(ws_test,density_test)

 #computing weights based on a Weibull ditribution for wind speed and a uniform distribution for ambient temperature 
 userweights = dweibull(testset[,1], shape = 2.25, scale = 6.5)*1 #assuming uniform distribution for ambient temperature
 userweights = userweights/sum(userweights) #normalizing the weights to ensure they sum to 1
 
 # Read the data
data1 = data1[1:100, ]
data2 = data2[1:100, ]
datalist = list(data1, data2)
xCol = c(2, 4) 
xCol.circ = NULL
yCol = 7
testCol = c(2, 4)
 
 # Compute power curve difference on the constructed testset
 output = ComparePCurve(data = datalist, xCol = xCol, yCol = yCol, testCol = testCol, testSet = testset) 
 
 #Compute weighted diff based on the calculated weights 
 weightedDiff = ComputeWeightedDifference(muDiff = output$muDiff, weights = userweights, base = output$mu1)

 #Compute statistically significant weighted diff
 weightedStatDiff = ComputeWeightedDifference(muDiff = output$muDiff, weights = userweights, base = output$mu1, statDiff = TRUE, confBand = output$band)
 ```

### 3. funGP
 The function can be used to perform function comparison using Gaussian process and hypothesis testing

 *Function :*

 *funGP (datalist, xCol, yCol, confLevel = 0.95, testset, limitMemory = TRUE, opt_method = 'L-BFGS-B', sampleSize = list(optimSize = 500, bandSize = 5000), rngSeed = 1 )*

 ```R
 # Preparing the arguments
datalist = list(data1[1:100,], data2[1:100, ])
xCol = c(2, 4)
yCol = 7
confLevel = 0.95
testset = matrix(c(7.2, 1.14, 12.3, 1.16), nrow = 2, ncol = 2, byrow = TRUE)
limitMemory = TRUE
opt_method = 'L-BFGS-B'
sampleSize = list(optimSize = 500, bandSize = 5000) 
rngSeed = 1 

 # Executing the function
function_diff = funGP(datalist, xCol, yCol, confLevel, testset, limitMemory, opt_method, sampleSize, rngSeed)
 ```
 
### 4. tempGP
  A Gaussian process based power curve model which explicitly models the temporal aspect of the power curve. The model consists of two parts: f(x) and g(t).
  
  *Function :*

 *tempGP (trainX, trainY, trainT = NULL)*
 
 ```R
 data = DSWE::data1
 trainindex = 1:5000 #using the first 5000 data points to train the model
 traindata = data[trainindex,]
 xCol = c(2:6) #input variable columns
 yCol = 7 #response column
 tCol = 1 #column with the time indices
 trainX = as.matrix(traindata[,xCol])
 trainY = as.numeric(traindata[,yCol])
 trainT = as.numeric(traindata[,tCol])
 
 ## Two ways to call the tempGP function
 # 1. Using user defined trainT
 tempGPObject = tempGP(trainX, trainY, trainT)
 
 # 2. Generate time indices internally as the sequence of integers starting from 1.
 tempGPObject = tempGP(trainX, trainY)
 
 ## Use the predict() method to get predictions from the learned model.
 testdata = data[5001:10000,] # defining test data as the next 5000 data points after train indices
 
 ## Predict only the function f(x) and ignore temporal component g(t)
 testX = as.matrix(testdata[,xCol,drop = F])
 testY = as.numeric(testdata[,yCol])
 predF = predict(tempGPObject, testX)
 rmseF = sqrt(mean((testY - predF)^2)) #rmse 
 cat('RMSE using just f(x):',round(rmseF,3),'\n')
 
 ## Or predict both f(x) and g(t) using a rolling window data update with the help of updateData() method for tempGP objects. 
 
 predY = rep(0,nrow(testdata)) #vector to store the rolling predictions
 for (i in 1:nrow(testdata)){
  testX_i = as.matrix(testdata[i,xCol,drop = F]) #input variables for time point i; replace with forecast when actual data not available.
  testT_i = testdata[i,tCol] #time index for i
  predY[i] = predict(tempGPObject, testX_i, testT_i) #predict both f(x) and g(t)
  
  #After time point i, the data for time point i would be available. Update the data and residuals in the tempGP object.
  
  tempGPObject = updateData(tempGPObject, newX = testX_i, newY = testY[i], newT = testT_i)
 }
 rmseY = sqrt(mean((testY - predY)^2)) #rmse 
 cat('RMSE using rolling update and f(x) + g(t):',round(rmseY,3),'\n')
 ```
 
 

### 5. KnnPCFit
The function can be used to model the data using user supplied arguments, a knn model is returned as an end result. It can also be used to get the best feature subset, if subsetSelection is set TRUE

*Function :*

*KnnPCFit(data, xCol, yCol, subsetSelection = FALSE)*

```R
# Preparing the arguments
data = data1
xCol = c(2, 4)
yCol = 7
subsetSelection = FALSE

# Executing the function
knn_model = KnnPCFit(data, xCol, yCol, subsetSelection)
```
### 6. KnnPredict
The function can be used to evaluate a prediction on a new test point using model generated using KnnPCFit

*Function :*

*KnnPredict(knnMdl, testData)*

```R
# Preparing the arguments
data = data1
xCol = c(2, 4)
yCol = 7
subsetSelection = FALSE

knn_model = KnnPCFit(data, xCol, yCol, subsetSelection)
testData = data[1:100, ]

# Executing the function
prediction = KnnPredict(knn_model, testData)
```

### 7. KnnUpdate
The function can be used to update the knn model whenever new data in available

*Function :*

*KnnUpdate = function(knnMdl, newData)*

```R
# Preparing the arguments
data = data1
xCol = c(2, 4)
yCol = 7
subsetSelection = FALSE

knn_model = KnnPCFit(data, xCol, yCol, subsetSelection)
newData = data[500:1000, ]

# Executing the function
knn_newmodel = KnnUpdate(knn_model, newData)
```

### 8. AMK
The function can be used to model the data by using user supplied arguments. It uses a kernel to assign weights to every training data points, the bandwidth of kernel (bw) can be provided as vector of values or character 'dpi' and 'dpi-gap'. If provided character input, the bandwidths are computed internally

*Function :*

*AMK(trainX, trainY, testX, bw = 'dpi_gap', nMultiCov = 3, fixedCov = c(1, 2), cirCov = NA )*

```R
# Preparing the arguments
data = data1
trainX = data[, c(2, 4)]
trainY = data[, 7]
testX = data[100:140, c(2, 4)]
bw = 'dpi_gap'
nMultiCov = 2
fixedCov = NULL
cirCov = NA

# Executing the function
AMK_prediction = AMK(trainX, trainY, testX, bw, nMultiCov, fixedCov, cirCov)
```
**Note :-** In case an error such as a non-finite bandwidth is generated upon adding a covariate, please remove such covariates

### 9. BayesTreePCFit
The function can be used to model the data by using user supplied arguments. It uses tree based method to model the supplied data set

*Function :*

*BayesTreePCFit(trainX, trainY, testX, nTree = 50)*

```R
# Preparing the arguments
data = data1
trainX = data[, c(2, 4)]
trainY = data[, 7]
testX = data[100:110, c(2, 4)]

# Executing the function
Bart_prediction = BayesTreePCFit(trainX, trainY, testX)
```
### 10. SplinePCFit
The function can be used to model the data by using user supplied arguments. It uses spline based method to model the data set. Users can leverage the modelFormula argument to model data set using interactions among features

*Function :*

*SplinePCFit(data, xCol, yCol, testP, modelFormula = NULL)*

```R
# Preparing the arguments
data = data1
xCol = c(2, 4)
yCol = 7
testX = data[100:110, ]

# Executing the function
Spline_prediction = SplinePCFit(data, xCol, yCol, testX)
```
### 11. SvmPCFit
The function can be used to model the data by using user supplied arguments. It uses support vector machine to model the data set.

*Function :*

*SvmPCFit(trainX, trainY, testX, kernel = 'radial')*

```R
# Preparing the arguments
data = data1
trainX = data[, c(2, 4)]
trainY = data[, 7]
testX = data[100:110, c(2, 4)]

# Executing the function
Svm_prediction = SvmPCFit(trainX, trainY, testX)
```
### 12. CovMatch
The function can be used to match different data sets. It can only be used to match two different data set at one time. If priority argument is set to FALSE, which is default, the feature columns provided are used in the same order in matching, else computes the covariates matching sequence

*Function :*

*CovMatch(data, xCol, xCol.circ = NULL, thrs = 0.2, priority = FALSE)*

```R
# Preparing the arguments
data1 = data1[1:100, ]
data2 = data2[1:100, ]

data = list(data1, data2)
xCol = c(2, 4)
xCol.circ = NULL
thrs = c(0.1, 0.1)
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
