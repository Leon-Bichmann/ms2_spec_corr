library(protViz)
library("mzR") ## the software package
library("msdata") ## the data package

#####################
##### Variables #####
sequence<-'SRFLSQLDK'
n<-	1
label_pos<-5
label_weigth<-1
ftol<-0.02


######################
#####  Functions #####
cos.sim <- function(A,B)
{
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}

######################
#####    Main    #####

#read measured spectra
file <- '/Users/bichmann/Projects/synpeptides/SRFLSQLDK_Mel15.mzML'
mz <- openMSfile(file)
h<-header(mz)
spec_nr<-which((h$acquisitionNum==n) & (h$msLevel==2))
p<-as.data.frame(spectra(mz,scans=spec_nr))
names(p)<-c('mZ','intensity')
plot<-peakplot(sequence,spec=p, ion.axes = F, itol = ftol)

# add charge 2 fragments
fs<-plot$fragmentIon
colnames(fs)<-c("b+","y+","c+","z+")
plot$fragmentIon$c<-(plot$fragmentIon$c+1)/2
plot$fragmentIon$z<-(plot$fragmentIon$z+1)/2
plot$fragmentIon$b<-(plot$fragmentIon$b+1)/2
plot$fragmentIon$y<-(plot$fragmentIon$y+1)/2
tfs<-cbind(plot$fragmentIon, fs)
plot<-peakplot(sequence,spec=p, ion.axes = F, fi = tfs, itol = ftol)

#isotope shift b and c ions
plot$fragmentIon["b"][label_pos:nchar(sequence),]<-plot$fragmentIon["b"][label_pos:nchar(sequence),]+label_weigth
plot$fragmentIon["c"][label_pos:nchar(sequence),]<-plot$fragmentIon["c"][label_pos:nchar(sequence),]+label_weigth
plot$fragmentIon["b+"][label_pos:nchar(sequence),]<-plot$fragmentIon["b+"][label_pos:nchar(sequence),]+label_weigth
plot$fragmentIon["c+"][label_pos:nchar(sequence),]<-plot$fragmentIon["c+"][label_pos:nchar(sequence),]+label_weigth

#isotope shift y and z ions
plot$fragmentIon["y"][(nchar(sequence)-label_pos):nchar(sequence),]<-plot$fragmentIon["y"][(nchar(sequence)-label_pos):nchar(sequence),]+label_weigth
plot$fragmentIon["z"][(nchar(sequence)-label_pos):nchar(sequence),]<-plot$fragmentIon["z"][(nchar(sequence)-label_pos):nchar(sequence),]+label_weigth
plot$fragmentIon["y+"][(nchar(sequence)-label_pos):nchar(sequence),]<-plot$fragmentIon["y+"][(nchar(sequence)-label_pos):nchar(sequence),]+label_weigth
plot$fragmentIon["z+"][(nchar(sequence)-label_pos):nchar(sequence),]<-plot$fragmentIon["z+"][(nchar(sequence)-label_pos):nchar(sequence),]+label_weigth

#store intensitys of respective fragments
tfs<-plot$fragmentIon
plot<-peakplot(sequence,spec=p, ion.axes = F, fi = tfs, itol = ftol)
Ints_measured<-p$intensity[plot$idx]

#########################
#read synthesized spectra
file <- '/Users/bichmann/Projects/synpeptides/SRFLSQLDK_syn_HCD.mzML'
mz <- openMSfile(file)
h<-header(mz)
spec_nr<-which((h$acquisitionNum==n) & (h$msLevel==2))
p<-as.data.frame(spectra(mz,scans=spec_nr))
names(p)<-c('mZ','intensity')
plot<-peakplot(sequence,spec=p, ion.axes = F, itol = ftol)

# add charge 2 fragments
fs<-plot$fragmentIon
colnames(fs)<-c("b+","y+","c+","z+")
plot$fragmentIon$c<-(plot$fragmentIon$c+1)/2
plot$fragmentIon$z<-(plot$fragmentIon$z+1)/2
plot$fragmentIon$b<-(plot$fragmentIon$b+1)/2
plot$fragmentIon$y<-(plot$fragmentIon$y+1)/2
tfs<-cbind(plot$fragmentIon, fs)

#store intensitys of respective fragments
plot<-peakplot(sequence,spec=p, ion.axes = F, fi = tfs, itol = ftol)
Ints_synthesized<-p$intensity[plot$idx]

#compute cosine similarity
cos.sim(Ints_measured,Ints_synthesized)
