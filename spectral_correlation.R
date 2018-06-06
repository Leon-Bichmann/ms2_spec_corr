library(protViz)
library("mzR") ## the software package
library("msdata") ## the data package

#####################
##### Variables #####

sequence<-'SYFPEITHI' # Peptide Sequence of Interest
n<-	1                 # Spectrum Number (=1 if extracted spectra)
label_pos<-2          # Labelled residue index eg. Y in SYFPEITHI example
label_weigth<-1       # Mass delta of isotope label eg. +1
ftol<-0.5             # fragment tolerance in Dalton for assigning y, b, z and c ions


######################
#####  Functions #####

cos.sim <- function(A,B)
{
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}

######################
#####    Main    #####

#read measured spectra
file <- '/Users/bichmann/Projects/synpeptides/SYFPEITHI.mzML'
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
file <- '/Users/bichmann/Projects/synpeptides/SYFPEITHI_syn_HCD.mzML'
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
