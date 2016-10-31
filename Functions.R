######################
#
#   Grep name in x
#
#####################
specgrep <- function(x,name) {
  x=x[,grep(name, names(x))]
  return(x)
}


######################
#
#   Remove name in x
#
#####################
removegrep <- function(x,name) {
  x=x[,-grep(name, names(x))]
  return(x)
}


######################
#
#   BlankFilter - find features which are "highly" present in blank (based on a ratio cutoff)
#   to.remove=AdvancedBlankFilter(Blanks,samples,0.01)
#
#####################
BlankFilter <- function(blanks, samples, cutoff) {
  blanks[is.na(blanks)] <- 0
  samples[is.na(samples)] <- 0
  
  blanks <- apply(blanks,1,median,na.rm=TRUE)
  samples <- apply(samples,1,max,na.rm=TRUE)
  
  to.remove <- which(blanks/samples >= cutoff)
  return(to.remove)
}


######################
#
#   remove.batchfeatures - find "batch specific" features 
#   to.remove=remove.batchfeatures(samples)
#
#####################
remove.batchfeatures <- function(samples) {
  B1_samples=samples[,grep("B1", names(samples))]
  B2_samples=samples[,grep("B2", names(samples))]
  B3_samples=samples[,grep("B3", names(samples))]
  B4_samples=samples[,grep("B4", names(samples))]
  
  idx_1=coverage(B1_samples,0.8)
  idx_2=coverage(B2_samples,0.8)
  idx_3=coverage(B3_samples,0.8)
  idx_4=coverage(B4_samples,0.8)
  
  presence=matrix(0,dim(samples)[1],4)
  presence[idx_1,1] = 1
  presence[idx_2,2] = 1
  presence[idx_3,3] = 1
  presence[idx_4,4] = 1
  idx=0
  for (i in 1:dim(samples)[1]) {
    row=presence[i,]
    if (length(row[row==1])==1) {
      idx = c(idx,i)
    }
  }
  
  idx=idx[idx>0]
  return(idx)
}


######################
#
#   coverage - compute features with coverage of cutoff or above
#   features=coverage(samples,1)
#
#####################
coverage <- function(intensities, coverage) {
  amount = dim(intensities)
  cutoff = amount[2]*coverage
  features=rep(0,amount[1])
  
  for (i in 1:amount[1]) {
    temp=intensities[i,]
    a=!is.na(temp)
    a=which(a)
    if (length(a)>=cutoff) {
      features[i]=i
    }
  }
  features=features[which(features>0)]
  return(features)
}


######################
#
#   extract.names - extract unique names
#   names=extract.names(samples)
#
#####################
extract.names <- function(x){
  names=gsub("intensity_","",names(x))
  names=gsub("_Rep1","",names)
  names=gsub("_Rep2","",names)
  names=gsub("_Rep3","",names)
  
  names=gsub("_B1","",names)
  names=gsub("_B2","",names)
  names=gsub("_B3","",names)
  names=gsub("_B4","",names)
  names=gsub("_2","",names)
  names=unique(names)
  return(names)
}


######################
#
#   remove.narows - removes rows/features which are completely missing (NA)
#   samples=remove.narows(samples)
#
#####################
remove.narows<-function(X)  {
  ind <- apply(X, 1, function(x) all(is.na(x)))
  return(ind)
}


######################
#
#   Outersect
#
#####################
outersect <- function(x, y) {
  sort(c(x[!x%in%y],
         y[!y%in%x]))
}

