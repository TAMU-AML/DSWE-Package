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
                           sample=sample, size=size)
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