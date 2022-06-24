syncSize = function(data, powercol, timecol=0, xcol, method="minimum power"){
 
   # Checks whether the provided data set is a list or not
  
  if(!is.list(data)){
    stop('data must be a list containing pair of data sets')
  }
  if(length(data) != 2){
    stop('The number of data sets to synchronize should be equal to two')
  }
  
  if(!is.numeric(xcol)) {
    stop('The xcol must be provided as a numeric/vector')
  }
  
  #powercol is only used when method set to 'minimum power'
    if (method == 'minimum power'){
      if(!is.numeric(powercol)){
        stop('The powercol must be provided as a numeric')
      }else if(length(powercol) > 1){
        stop('The powercol should be a single numeric value')
      }
      if(powercol==0){
        stop('The powercol cannot be zero')
      }
    }
    
    if(method=="time" & timecol==0) {
      stop('Provide value for timecol and make sure data sets have the same time format')
    }
    
    if (method == "time" && nchar(data[[1]][1,timecol]) != nchar(data[[2]][1,timecol])) {
      stop('Make sure data sets have same the time format')
    }
    
    ## Quick clean and prepare the data
    library(dplyr)
    Data1 = data[[1]]; Data2 = data[[2]]
    if (timecol > 0){
      names(Data1)[timecol]=names(Data2)[timecol]="Timestamp" 
    }
    names(Data1)[powercol]=names(Data2)[powercol]="Power"
    col=c(timecol,powercol,xcol)
    col=unique(col)
    Data1 = Data1[,col]; Data2 = Data2[,col]
    Data1 = Data1[complete.cases(Data1),]; Data2 = Data2[complete.cases(Data2),]
    Data2.sync = Data2; Data1.sync = Data1
    n1 = nrow(Data1); n2=nrow(Data2)
    del = abs(n1-n2)
    if (del>0){
      ## Synchronize with "minimum power"
      if (method == "minimum power"){
        if (n1 > n2){
          idx = which(Data1$Power<=0)
          if (length(idx) > del){Data1.sync=Data1[-idx[1:del],]} else {
            for (m in length(idx):del){
              i = length(idx)+1
              a = min(Data1$Power[-idx])
              b = del - length(idx)
              new = which(Data1$Power == a)
              if (length(new) >= b) {idx[i:del]=new[1:b]} else {
                idx[i:(i+length(new)-1)]=new }
              if (length(idx)==del){break}
            }
            Data1.sync=Data1[-idx,]}
        }
        if (n2>n1){
          idx = which(Data2$Power<=0)
          if (length(idx) > del){Data2.sync=Data2[-idx[1:del],]} else {
            for (m in length(idx):del){
              i = length(idx)+1
              a = min(Data2$Power[-idx])
              b = del - length(idx)
              new = which(Data2$Power == a)
              if (length(new) >= b) {idx[i:del]=new[1:b]} else {
                idx[i:(i+length(new)-1)]=new }
              if (length(idx)==del){break}
            }
            Data2.sync=Data2[-idx,]}
        }
      }
      ## Synchronize with "time"
      if (method == "time"){
        idx1 = unique(Data1$Timestamp)
        idx2 = unique(Data2$Timestamp)
        timestamp = c(idx1,idx2)
        timestamp = timestamp[duplicated(timestamp)]
        Data1.sync = Data1 %>% filter(Timestamp %in% timestamp)
        Data2.sync = Data2 %>% filter(Timestamp %in% timestamp)
        for (z in 1:20){
          idx1 = which(duplicated(Data1.sync$Timestamp)==TRUE)
          idx2 = which(duplicated(Data2.sync$Timestamp)==TRUE)
          if (length(idx1)==0 & length(idx2)==0){break} else{
            Data1.sync = Data1.sync[-idx1,]
            Data2.sync = Data2.sync[-idx2,]
          }
        }
      }
      ## Synchronize with "random"
      if(method == "random"){
        if(n1>n2){
          idx=sample(1:nrow(Data1),del)
          Data1.sync = Data1[-idx,]} else {if(n2>n1){
            idx=sample(1:nrow(Data2),del)
            Data2.sync = Data2[-idx,]}
          }
      } 
    }
    Data = list(Data1=Data1.sync,Data2=Data2.sync)
    return(Data) 
}