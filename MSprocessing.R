

setwd("C:\\path\\to\\working\\directory")
source('Functions.R')

### Initial look at the data
library(rgl)

plotdata=log2(samples)
plotdata[is.na(plotdata)]=0
pca <- prcomp(t(plotdata))
name=names(samples)
name=gsub("intensity_","",name)

plot3d(pca$x[,c(1,2,3)], pch=20,size=10)

### Blankfiltering - removing contaminants
blanks=specgrep(samples,"BLANK")
run <- run[-grep("BLANK",run)]
samples=removegrep(samples,"BLANK")

to.remove=BlankFilter(blanks,samples,0.01)
samples=samples[-to.remove,]
cat("Removed by BlankFilter: ",length(to.remove), "\n")
rm(to.remove)

### Remove batch specific features
to.remove=remove.batchfeatures(samples)
if(length(to.remove)>0) {samples=samples[-to.remove,]
cat("Batch specific removed: ",length(to.remove), "\n")} else {cat("No batch specific features")}
rm(to.remove)

### log2 transformation
samples = log2(samples)

### Outlier exclusion
TIC<-colSums(samples, na.rm=T)
par(mfrow=c(1,2))
boxplot(TIC, ylab="Extracted peaks area")
hist(TIC, breaks=100, main="", xlab="Extracted peaks area")
samples <- samples[,!(colnames(samples) %in% c(colnames(samples)[TIC<mean(TIC)*0.7]))]

### Outlier exclusion by replicate correlation
samples_temp = samples
drugs = extract.names(samples_temp)
cutoff = 0.80

y = specgrep(samples_temp,drugs[1])
c = cor(y, use="complete.obs") 
c[lower.tri(c, diag=TRUE)]=0
ind = which(c>cutoff,arr.ind=TRUE)
ind = unique(c(ind))
samples = y[,ind]

excluded = outersect(seq(1:dim(y)[2]),ind)
names=names(y)
if (length(excluded>0)){
  cat("The following were excluded")
  cat(names[excluded], "\n", "Min cor: ",min(c[c>0]),"Max cor: ",max(c),"\n")}

for (i in 2:length(drugs)) {
  y = specgrep(samples_temp,drugs[i])
  c = cor(y, use="complete.obs") 
  c[lower.tri(c, diag=TRUE)]=0
  c = round(c,2)
  print(c[c>0])
  ind = which(c>=cutoff,arr.ind=TRUE)
  ind = unique(c(ind))
  samples = cbind(samples,y[,ind])
  excluded = outersect(seq(1:dim(y)[2]),ind)
  names=names(y)
  
  if (length(excluded>0)){cat(names[excluded], "\n", "Min cor: ",min(c[c>0]),"Max cor: ",max(c),"\n")}
}

### Remove NA rows
ind = remove.narows(samples)
if (any(ind)) {  samples <- samples[!ind, ]}

### Normalize by LOESS
library(limma)
par(mfrow=c(2,1))
boxplot(samples)
samples<-data.frame(normalizeCyclicLoess(samples))
boxplot(samples)

### Save NA positions for later reintroduction
samples[samples==0]=NA
NAs <- which(is.na(samples), arr.ind=TRUE)

### Sort in runorder
run = run[-which(run %in% outersect(run,names(samples)))]
samples = samples[,run]
cat("Cut happens at position: ", grep("^Mebendazole_B4_Rep2$",names(samples)))

### Write to file
samples[is.na(samples)] = 0
write.table(samples,"preprocessed_data.xls",sep = "\t", row.names = F)


########################################################################################################################
#                                                                                                                      #
#                                     Perform OOS-DA to remove unwanted variation                                      #
#                                                                                                                      #
########################################################################################################################

### Reintroduce NAs
### corrected = OOS-DA output
temp = samples[,which(names(samples) %in% names(corrected))]
temp[temp==0]=NA
NAss <- which(is.na(temp), arr.ind=TRUE)

corrected = corrected[,names(temp)]
corrected[NAss] = NA
rm(temp)

### Remove NA rows
ind = remove.narows(corrected)
if (any(ind)) {  corrected <- corrected[!ind, ]}

### Merge replicates
temp = corrected 
rm(drugs)
drugs = extract.names(temp)
drugs = drugs[-grep("Contro",drugs)]
merged = specgrep(temp,"Control") # exclude controls from replicate merge
no = ncol(merged)

for (i in 1:length(drugs)) {
  y = specgrep(temp,drugs[i])
  y = apply(y, 1, median, na.rm=TRUE) 
  merged = cbind(merged,y)
}
names = c(names(merged[1:no]),drugs)
names(merged)=names

### Sort by mechanism of action
C <- grep("Control", names(merged))
T <- outersect(seq(1:ncol(merged)), C)
merged = merged[,c(C,T)]
cat("Controls: ", length(C), "Treated: ", length(T))

merged[is.na(merged)] = 0
write.table(merged,"final_data.xls",sep = "\t", row.names = F)
