# MIT License
# 
# Copyright (c) 2022 Effi Latiffianti and Yu Ding
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#   
#   The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

#' @title Power imputation
#' 
#' @description Good power curve modeling requires valid power values in the region between cut-in and cut-out wind speed. However, when turbine is not operating, the power production will be recorded as zero or negative. This function replaces those values with predicted values obtained from the estimated tempGP power curve model using one input variable - the wind speed. 
#'
#' @param data A list of two data sets that require imputation.
#' @param powercol A numeric stating the column number of power production.
#' @param vcol  A numeric stating the column number of wind speed.
#' @param vrange A vector of cut-in, rated, and cut-out wind speed.
#' @param rated.power A numerical value stating the wind turbine rated power. 
#' @param sample A boolean (TRUE/FALSE) indicating whether to use sample or the whole data sets to train the power curve. 
#' @param size A numeric stating the size of sample when \code{sample = TRUE}. Default value is 2500. It is only used when \code{sample = TRUE}.
#' 
#' @return a list containing datasets with the imputed power. 
#' 
#' @importFrom stats complete.cases
#'  
#' @examples 
#' 
#' data = list(data1[1:100,], data2[1:120, ])
#' powercol = 7
#' vcol = 2
#' vrange = c(5,12,25)
#' rated.power = 100
#' sample = FALSE
#' 
#' imputed.dat = imptPower(data, powercol, vcol, vrange, rated.power, sample)
#' 
#' @references Latiffianti, E, Ding, Y, Sheng, S, Williams, L, Morshedizadeh, M, Rodgers, M (2022). "Analysis of leading edge protection application on wind turbine performance through energy and power decomposition approaches". Wind Energy. 2022; 1-19. \doi{10.1002/we.2722}.
#' 
#' @export

imptPower = function(data, powercol, vcol, vrange, rated.power=NULL, sample = TRUE, size = 2500){
  
  if(!is.list(data)){
    stop('data must be a list containing pair of data sets')
  }
  
  if(length(data) != 2){
    stop('The number of data sets for imputation should be equal to two')
  }
  
  if(!is.numeric(powercol)){
    stop('The powercol must be provided as a numeric')
  }else if(length(powercol) > 1){
    stop('The powercol should be a single numeric value')
  }
  
  if(!is.numeric(vcol)){
    stop('The vcol must be provided as a numeric')
  }else if(length(vcol) > 1){
    stop('The vcol should be a single numeric value')
  }
  
  if(!is.numeric(vrange)){
    stop('The vrange must be provided as a vector contains three numeric values')
  }else if(length(vrange) != 3){
    stop('The vrange should provide three numeric values of cut-in, rated, and cut-out wind speed')
  }
  
  if (!inherits(sample, "logical")){
    stop('sample should either be TRUE or FALSE')
  }
  
  if(!is.numeric(size)){
    stop('The size must be provided as a numeric. It is the number sample for power curve modeling.')
  }else if(length(size) > 1){
    stop('The size should be a single numeric value')
  }
  
  Data1 = data[[1]][,c(powercol,vcol)]
  Data2 = data[[2]][,c(powercol,vcol)]
  names(Data1)=names(Data2)=c("Power","Windspeed")
  Data1$obs = c(1:nrow(Data1)); Data2$obs = c(1:nrow(Data2))
  
  dat.max1 = Data1$Power[which(Data1$Windspeed > (vrange[2]) & Data1$Windspeed < vrange[3])]
  dat.max2 = Data2$Power[which(Data2$Windspeed > (vrange[2]) & Data2$Windspeed < vrange[3])]
  if(min(length(dat.max1),length(dat.max2))== 0 & is.null(rated.power)==TRUE) {
    stop('provide value for rated.power, no data above rated wind speed is available')
  }
  
  if (length(dat.max1)==0) {max1 = rated.power} else {max1 = mean(dat.max1)}
  if (length(dat.max2)==0) {max2 = rated.power} else {max2 = mean(dat.max2)}
  max.power = c(max1,max2)
  
  dat1 = Data1[which(Data1$Windspeed >= vrange[1] & Data1$Windspeed <= vrange[2]),]
  dat2 = Data2[which(Data2$Windspeed >= vrange[1] & Data2$Windspeed <= vrange[2]),]
  dat1 = dat1[which(dat1$Power > 0),]
  dat2 = dat2[which(dat2$Power > 0),]
  dat1 = dat1[complete.cases(dat1),]
  dat2 = dat2[complete.cases(dat2),]
  
  if (nrow(dat1) == 0 || nrow(dat2) == 0) { 
    stop('not enough data point for predicting the imputed values')}
  
  if(sample){
    idx1 = sample(c(1:nrow(dat1)),min(size,nrow(dat1)))
    idx2 = sample(c(1:nrow(dat2)),min(size,nrow(dat2)))
    if (nrow(dat1) < size){idx1 = 1:nrow(dat1)}
    if (nrow(dat2) < size){idx2 = 1:nrow(dat2)}
    PC1 = tempGP(dat1$Windspeed[idx1],dat1$Power[idx1])
    PC2 = tempGP(dat2$Windspeed[idx2],dat2$Power[idx2])
  } else {
    PC1 = tempGP(dat1$Windspeed,dat1$Power)
    PC2 = tempGP(dat2$Windspeed,dat2$Power)
  }
  
  dat.list = list(Data1[complete.cases(Data1$Windspeed),],Data2[complete.cases(Data1$Windspeed),])
  PC = list(PC1,PC2)
  for (i in 1:2) {
    m = which(dat.list[[i]]$Power<=0)
    dat.list[[i]]$Power[m]=0
    j1 = which(dat.list[[i]]$Windspeed[m] > vrange[1] & dat.list[[i]]$Windspeed[m]<= vrange[2])
    j2 = which(dat.list[[i]]$Windspeed[m] > vrange[2] & dat.list[[i]]$Windspeed[m]<= vrange[3])
    imput.j1 = predict(PC2,dat.list[[i]]$Windspeed[m[j1]])
    imput.j1[which(imput.j1>max.power[i])]=max.power[i]
    imput.j1[which(imput.j1<0)]=0
    dat.list[[i]]$Power[m[j1]]= imput.j1
    dat.list[[i]]$Power[m[j2]]= max.power[i]
  }
  Data1 = data[[1]][dat.list[[1]]$obs,]
  Data2 = data[[2]][dat.list[[2]]$obs,]
  Data1[,powercol] = dat.list[[1]]$Power
  Data2[,powercol] = dat.list[[2]]$Power
  Data = list(Data1 = Data1, Data2=Data2)
  return(Data)
}
