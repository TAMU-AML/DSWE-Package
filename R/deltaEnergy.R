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

#' @title Energy decomposition for wind turbine performance comparison
#' 
#' @description Energy decomposition compares energy production from two datasets and separates it into turbine effects (deltaE.turb) and weather/environment effects (deltaE.weather). 
#'
#' @param data A list of two data sets to be compared. A difference  is always computed as (data2 - data1).
#' @param powercol A numeric stating the column number of power production.
#' @param timecol A numeric stating the column number of data time stamp. Default value is zero. A value other than zero should be provided when \code{sync.method = 'time'}. 
#' @param xcol A numeric or vector stating the column number(s) of power curve input covariates/features (environmental or weather variables are recommended).
#' @param sync.method A string specifying data synchronization method. Default value \code{'minimum power'}; other options include \code{'time'} and \code{'random'}.
#' @param imput A boolean (TRUE/FALSE) indicating whether power imputation should be performed before calculating energy decomposition. The recommended and default value is TRUE. Change to FALSE when data have been preprocessed or imputed before.#' @param vcol A numeric stating the column number of wind speed. It is required when \code{imput = TRUE}.
#' @param vcol A numeric stating the column number of wind speed.
#' @param vrange A vector of cut-in, rated, and cut-out wind speed. Values should be provided when \code{imput = TRUE}.
#' @param rated.power A numerical value stating the wind turbine rated power. 
#' @param sample A boolean (TRUE/FALSE) indicating whether to use sample or the whole data sets to train the power curve to be used for power imputation. Default value is TRUE. It is only used when \code{imput = TRUE}.
#' @param size A numeric stating the size of sample when \code{sample = TRUE}. Default value is 2500. It is only used when \code{imput = TRUE} and \code{sample = TRUE}.
#' @param timestamp.min A numerical value stating the resolution of the datasets in minutes. It is the difference between two consecutive time stamps at which data were recorded. Default value is 10. 
#' 
#' @return a list containing :
#'  \itemize{
#'   \item deltaE.turb - A numeric, % of energy difference representing the turbine effects. 
#'   \item deltaE.weather - A numeric, % of energy difference representing the weather effects.
#'   \item deltaE.hat - A numeric, % of total energy difference estimated.
#'   \item deltaE.obs - A numeric, % of total energy difference observed.
#'   \item estimated.energy - A numeric vector of the total energy calculated from each of f1(x2), f1(x1),  f2(x2), f1(x2). If power is in kW, these values will be in kWh.
#'   \item data - A list of two datasets used to calculate energy decomposition, i.e. synchronized. When \code{imput = TRUE}, the power column is the result from imputation.
#' }
#'  
#' @examples 
#' 
#' data = list(data1[1:210,], data2[1:400, ])
#' powercol = 7
#' timecol = 1
#' xcol = c(2:6)
#' sync.method = 'time'
#' imput = TRUE
#' vcol = 2
#' vrange = c(5,12,25)
#' rated.power = 100
#' sample = FALSE
#' Decomposition = deltaEnergy(data, powercol, timecol, xcol, sync.method, imput,
#' vcol, vrange, rated.power, sample)
#' 
#' @references Latiffianti, E, Ding, Y, Sheng, S, Williams, L, Morshedizadeh, M, Rodgers, M (2022). "Analysis of leading edge protection application on wind turbine performance through energy and power decomposition approaches". Wind Energy. 2022; 1-19. doi:10.1002/we.2722.<\url{https://onlinelibrary.wiley.com/doi/10.1002/we.2722}>. 
#' 
#'@export

deltaEnergy = function(data, powercol, timecol=0, xcol, sync.method ="minimum power", imput=TRUE,
                vcol=NULL, vrange=NULL, rated.power = NULL, sample = TRUE, size = 2500, 
                timestamp.min = 10){
  
  if(!is.list(data)){
    stop('data must be a list containing pair of data sets')
  }
  if(length(data) != 2){
    stop('The number of data sets to synchronize should be equal to two')
  }
  
  if(sync.method=="time" & timecol==0) {
    stop('Provide value for timecol and make sure data sets have the same time format')
  }
  
  if (sync.method == "time" & nchar(data[[1]][1,timecol]) != nchar(data[[2]][1,timecol])) {
    stop('Make sure data sets have same the time format')}
  
  if (imput){
    if(!is.numeric(vcol)){
      stop('The vcol must be provided as a numeric')
    } else if (length(vcol) != 1 | vcol < 1 | vcol > ncol(data[[1]])){
      stop('The vcol must be a single numeric value within the range of data column number')
    }
  }
  
  
  if (imput) {
    if(!is.numeric(vrange)){
      stop('The vrange must be provided as a vector contains three numeric values')
    }else if(length(vrange) != 3){
      stop('The vrange should provide three numeric values of cut-in, rated, and cut-out wind speed')
    }
  }
  
  if (class(sample)!="logical"){
    stop('sample should either be TRUE or FALSE')
  }
  
  if(!is.numeric(size)){
    stop('The size must be provided as a numeric. It is the number sample for power curve modeling.')
  }else if(length(size) > 1){
    stop('The size should be a single numeric value')
  }
  
  if(!is.numeric(timestamp.min)){
    stop('timestamp.min must be provided as a numeric (in minute)')
  }else if(length(timestamp.min) > 1){
    stop('timestamp.min should be a single numeric value')
  }
  
  ## Synchronization and imputation
  xcol = unique(c(vcol,xcol))
  data.sync = syncSize(data=data, powercol=powercol, timecol=timecol, xcol=xcol, 
                       method=sync.method)
  if (timecol==0){
    powercol=1
    xcol=2:(1+length(xcol))
    vcol = xcol[1]
  } else {
    powercol = 2
    xcol = 3:(2+length(xcol))
    vcol = xcol[1]
  }
  
  if(imput){
    data.imput = imptPower(data=data.sync, powercol=powercol, vcol=vcol, vrange=vrange,
                           rated.power = rated.power, sample=sample, size=size)
    data.sync = data.imput
  }
  
  ## Data
  Data1 = data.sync[[1]]; Data2 = data.sync[[2]]
  Data = list(Data1 = Data1, Data2 = Data2)
  
  ## Decomposition
  ## Fit the power curves
    PC1 = tempGP(as.matrix(Data1[,xcol]),as.numeric(Data1$Power))
    PC2 = tempGP(as.matrix(Data2[,xcol]),as.numeric(Data2$Power))
  
  ## Estimate the energy
  F1X1 = sum(predict(PC1,Data1[,xcol])) * (timestamp.min/60)
  F1X2 = sum(predict(PC1,Data2[,xcol])) * (timestamp.min/60)
  F2X1 = sum(predict(PC2,Data1[,xcol])) * (timestamp.min/60)
  F2X2 = sum(predict(PC2,Data2[,xcol])) * (timestamp.min/60)
  
  ## Calculate Delta Energy
  E.turb1 = F2X1 - F1X1
  E.turb2 = F2X2 - F1X2
  E.obs1 = sum(Data1$Power) * (timestamp.min/60)
  E.obs2 = sum(Data2$Power) * (timestamp.min/60) 
  E.weather1 = F1X2 - F1X1
  E.weather2 = F2X2 - F2X1
  E.obs.avg = (E.obs2 + E.obs1)/2
  deltaE.obs = (E.obs2 - E.obs1)/E.obs.avg*100
  deltaE.hat = (F2X2 - F1X1)/E.obs.avg*100
  deltaE.turb =((E.turb1 + E.turb2)/2)/E.obs.avg*100
  deltaE.weather = ((E.weather1 + E.weather2)/2)/E.obs.avg*100
  estimated.energy = rbind("f1(x1)"=F1X1,"f1(x2)"=F1X2,"f2(x1)"=F2X1,"f2(x2)"=F2X2)
  returnList = list(deltaE.turb = deltaE.turb, deltaE.weather = deltaE.weather,
                    deltaE.hat = deltaE.hat, deltaE.obs = deltaE.obs, 
                    estimated.energy = estimated.energy, Data=Data)
  return(returnList)
}