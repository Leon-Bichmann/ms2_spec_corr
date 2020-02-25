library(protViz)
library("mzR") ## the software package
library("msdata") ## the data package

#####################
##### Variables #####

sequence<-'ESVGLTTAR' # Peptide Sequence of Interest
n<-	1                 # Spectrum Number (=1 if extracted spectra)
label_pos<- 7         # Labelled residue index eg. Y in SYFPEITHI example
label_weigth<-5       # Mass delta of isotope label eg. +1
ftol<-0.02             # fragment tolerance in Dalton for assigning y, b, z and c ions

scan=4833
scan_2=3738
#scan_3=8718

######################
#####  Functions #####

cos.sim <- function(A,B)
{
  return( sum(A*B)/sqrt(sum(A^2)*sum(B^2)) )
}


mass_shift <- function(plot, sequence, label_pos, label_weigth)
{
  #mass shift b and c ions
  plot$fragmentIon["b+"][label_pos:nchar(sequence),]<-plot$fragmentIon["b+"][label_pos:nchar(sequence),]+label_weigth
  plot$fragmentIon["c+"][label_pos:nchar(sequence),]<-plot$fragmentIon["c+"][label_pos:nchar(sequence),]+label_weigth
  plot$fragmentIon["b"][label_pos:nchar(sequence),]<-plot$fragmentIon["b"][label_pos:nchar(sequence),]+label_weigth/2
  plot$fragmentIon["c"][label_pos:nchar(sequence),]<-plot$fragmentIon["c"][label_pos:nchar(sequence),]+label_weigth/2
  
  #mass shift y and z ions
  plot$fragmentIon["y+"][(nchar(sequence)-(label_pos-1)):nchar(sequence),]<-plot$fragmentIon["y+"][(nchar(sequence)-(label_pos-1)):nchar(sequence),]+label_weigth
  plot$fragmentIon["z+"][(nchar(sequence)-(label_pos-1)):nchar(sequence),]<-plot$fragmentIon["z+"][(nchar(sequence)-(label_pos-1)):nchar(sequence),]+label_weigth
  plot$fragmentIon["y"][(nchar(sequence)-(label_pos-1)):nchar(sequence),]<-plot$fragmentIon["y"][(nchar(sequence)-(label_pos-1)):nchar(sequence),]+label_weigth/2
  plot$fragmentIon["z"][(nchar(sequence)-(label_pos-1)):nchar(sequence),]<-plot$fragmentIon["z"][(nchar(sequence)-(label_pos-1)):nchar(sequence),]+label_weigth/2
  
  return(plot)
}

######################
#####    Main    #####

#read synthesized spectra
#file <- '/Users/bichmann/Projects/Melanoma_MZML/ETLKPGTCVKR_syn_HCD.mzML'
#file <-'/Users/bichmann/Downloads/170915_AM_BD-ZH09_Prostate_W_8%_DDA_#1_400-650mz_msms25.mzML'
file <-'/Users/bichmann/Downloads/200211_AM_42isotopelabeled_crypticSynpeps_AbQc_20fmol-ul_DDA#1_msms30.mzML'
mz <- openMSfile(file)
h<-header(mz)
spec_nr<-which((h$acquisitionNum==scan) & (h$msLevel==2))
p<-as.data.frame(spectra(mz,scans=spec_nr))
names(p)<-c('mZ','intensity')
plot<-peakplot(sequence, spec=p, ion.axes = F, itol = ftol)

# add charge 2 fragments
fs<-plot$fragmentIon
colnames(fs)<-c("b+","y+","c+","z+")
plot$fragmentIon$c<-(plot$fragmentIon$c+1)/2
plot$fragmentIon$z<-(plot$fragmentIon$z+1)/2
plot$fragmentIon$b<-(plot$fragmentIon$b+1)/2
plot$fragmentIon$y<-(plot$fragmentIon$y+1)/2
tfs<-cbind(plot$fragmentIon, fs)
plot<-peakplot(sequence,spec=p, ion.axes = F, fi = tfs, itol = ftol)

plot<-mass_shift(plot, sequence, label_pos, label_weigth)

#store intensitys of respective fragments
tfs<-plot$fragmentIon
plot<-peakplot(sequence,spec=p, ion.axes = F, fi = tfs, itol = ftol)
labels<-plot$label[which(plot$mZ.Da.error<ftol & plot$mZ.Da.error>-ftol)]

#########################
#read measured spectra
#file <- '/Users/bichmann/Projects/synpeptides/FVPPTAISHF_Mel15.mzML'
file <-'/Users/bichmann/Downloads/161115_AM_BD-ZH08_Colon_W_10%_DDA_#3_400-650mz_msms37.mzML'

mz <- openMSfile(file)
h<-header(mz)
spec_nr<-which((h$acquisitionNum==scan_2) & (h$msLevel==2))
p_exp<-as.data.frame(spectra(mz,scans=spec_nr))
names(p_exp)<-c('mZ','intensity')
plot_exp<-peakplot(sequence,spec=p_exp, ion.axes = F, itol = ftol)

# add charge 2 fragments
fs<-plot_exp$fragmentIon
colnames(fs)<-c("b+","y+","c+","z+")
plot_exp$fragmentIon$c<-(plot_exp$fragmentIon$c+1)/2
plot_exp$fragmentIon$z<-(plot_exp$fragmentIon$z+1)/2
plot_exp$fragmentIon$b<-(plot_exp$fragmentIon$b+1)/2
plot_exp$fragmentIon$y<-(plot_exp$fragmentIon$y+1)/2
tfs<-cbind(plot_exp$fragmentIon, fs)

#store intensitys of respective fragments
plot_exp<-peakplot(sequence,spec=p_exp, ion.axes = F, fi = tfs, itol = ftol)
labels_exp<-plot_exp$label[which(plot_exp$mZ.Da.error<ftol & plot_exp$mZ.Da.error>-ftol)]

#compute cosine similarity
common_labels<-union(labels, labels_exp)
Ints_measured<-p_exp$intensity[plot_exp$idx][which(plot_exp$label %in% common_labels)]
Ints_measured<-replace(Ints_measured,Ints_measured==-Inf,0)
Ints_synthesized<-p$intensity[plot$idx][which(plot$label %in% common_labels)]
Ints_synthesized<-replace(Ints_synthesized,Ints_synthesized==-Inf,0)

spectral_angle<-1-(2*acos(cos.sim(Ints_measured,Ints_synthesized))/pi)
