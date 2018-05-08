#Loading data and creating dataframe of indices. 

source("C:/Users/Swu2/Documents/DiversityIndices.R")
library(readr)
library(ape)
library(ade4)
library(seqinr)
library(RColorBrewer)
library(dplyr)
library(Biostrings)

#vaccine 3D7 strain
vax3D7 <- "DENANANSAVKNNNNEEPSDKHIKEYLNKIQNSLSTEWSPCSVTCGNGIQVRIKPGSANKPKDELDYANDIEKKICKMEKCSSVFNVVNSSIGLI"

#load TEP data
setwd("T:/vaccine/rtss_malaria_sieve/")
TEP.data <- read_tsv("qdata/sequences/TEP.tsv")

#subset for subjects age 5-17 months, separated into placebo and vaccine
clinical.data <- read_csv("adata/RTSSclinicalData.csv")
clinical.data <- subset(clinical.data, ageCateg=="[5-17] months")
clinical.placebo <- subset(clinical.data, vaccine==0)
clinical.vaccine <- subset(clinical.data, vaccine==1)

#TEP data for all participants age 5-17 months, as well as separated into vaccine and placebo groups
TEP.all <- TEP.data[TEP.data$subject %in% clinical.data$id,]
TEP.placebo <- TEP.data[TEP.data$subject %in% clinical.placebo$id,]
TEP.vaccine <- TEP.data[TEP.data$subject %in% clinical.vaccine$id,]

#number of distinct subjects in TEP data for participants age 5-17 months
subjects.TEP <- TEP.all %>% distinct(subject)
subjects.TEP <- subjects.TEP$subject

#number of distinct subjects in TEP placebo group for participants age 5-17 months
subjects.placebo <- TEP.placebo %>% distinct(subject)
subjects.placebo <- subjects.placebo$subject

#number of distinct subjects in TEP vaccine group for participants age 5-17 months
subjects.vaccine <- TEP.vaccine %>% distinct(subject)
subjects.vaccine <- subjects.vaccine$subject

#=======================TEP dataframe of indices=====================================
TEP.df <- data.frame(matrix(ncol = 43, nrow = length(subjects.TEP)))
colnames(TEP.df) <- c("Sh", "N.Sh", 
                      "GS",
                      "H(0)","H(1)","H(2)", "H(3)",
                      "Ren(0)","Ren(2)","HC(3)",
                      "PolySites","UniqMuts", "MutFreq", "CorSegSites",
                      "Rao(samp)",
                      "Mf(raw)", "Mf(K80)", "Mf(Trs)", "Mf(Trv)", "Mf(Syn)", "Mf(Nsyn)", "Mf(dN2ds)", 
                      "Pi(raw)", "Pi(K80)", "Pi(Trs)", "Pi(Trv)", "Pi(Syn)", "Pi(Nsyn)", "Pi(dN2ds)",
                      "Ratio(raw)", "Ratio(K80)", "Ratio(Trs)", "Ratio(Trv)", "Ratio(Syn)", "Ratio(Nsyn)", "Ratio(dN2ds)",
                      "FAD(raw)", "FAD(K80)", "FAD(Trs)", "FAD(Trv)", "FAD(Syn)", "FAD(Nsyn)", "FAD(dN2ds)")
for(i in 1:length(subjects.TEP)){
  sub <- subjects.TEP[i]
  data <- TEP.all %>% filter(subject==sub) %>% select(nuc_sequence, reads)
  seq <- as.vector(data$nuc_sequence)
  count <- as.vector(data$reads)
  len <- nchar(seq[1])
  dst <- DNA.dist(seq)
  TEP.df[i,1:43] <- c(Shannon.entropy(count), N.Shannon.entropy(count), 
                      GiniSimpson(count), 
                      qD(count,0), qD(count,1), qD(count,2), qD(count,3), 
                      renyi(count,0), renyi(count,2), HCq(count,3),
                      poly.sites(seq)$psites, poly.sites(seq)$nmuts, 
                      MutationFreq(SortByMutations(seq,count)$nr,SortByMutations(seq,count)$nm,len), correct.seg.sites(seq, count),
                      rao(count,dst), 
                      sq.diversity2(seq,count,1)$Mf[1], sq.diversity2(seq,count,1)$Mf[2], sq.diversity2(seq,count,1)$Mf[3], sq.diversity2(seq,count,1)$Mf[4], 
                      sq.diversity2(seq,count,1)$Mf[5], sq.diversity2(seq,count,1)$Mf[6], sq.diversity2(seq,count,1)$Mf[7],
                      sq.diversity2(seq,count,1)$Pi[1], sq.diversity2(seq,count,1)$Pi[2], sq.diversity2(seq,count,1)$Pi[3], sq.diversity2(seq,count,1)$Pi[4],
                      sq.diversity2(seq,count,1)$Pi[5], sq.diversity2(seq,count,1)$Pi[6], sq.diversity2(seq,count,1)$Pi[7],
                      sq.diversity2(seq,count,1)$Ratio[1], sq.diversity2(seq,count,1)$Ratio[2], sq.diversity2(seq,count,1)$Ratio[3], sq.diversity2(seq,count,1)$Ratio[4],
                      sq.diversity2(seq,count,1)$Ratio[5], sq.diversity2(seq,count,1)$Ratio[6], sq.diversity2(seq,count,1)$Ratio[7],
                      sq.diversity2(seq,count,1)$FAD[1], sq.diversity2(seq,count,1)$FAD[2], sq.diversity2(seq,count,1)$FAD[3], sq.diversity2(seq,count,1)$FAD[4],
                      sq.diversity2(seq,count,1)$FAD[5], sq.diversity2(seq,count,1)$FAD[6], sq.diversity2(seq,count,1)$FAD[7])
}
rownames(TEP.df) <- subjects.TEP
vaccine <- ifelse(rownames(TEP.df) %in% TEP.vaccine$subject, "vaccine", "placebo")
#save as csv file
write.csv(TEP.df, "C:/Users/Swu2/Documents/TEP.df")

#==================Placebo dataframe of indices================================
placebo.df <- data.frame(matrix(ncol = 43, nrow = length(subjects.placebo)))
colnames(placebo.df) <- c("Sh", "N.Sh", 
                      "GS",
                      "H(0)","H(1)","H(2)", "H(3)",
                      "Ren(0)","Ren(2)","HC(3)",
                      "PolySites","UniqMuts", "MutFreq", "CorSegSites",
                      "Rao(samp)",
                      "Mf(raw)", "Mf(K80)", "Mf(Trs)", "Mf(Trv)", "Mf(Syn)", "Mf(Nsyn)", "Mf(dN2ds)", 
                      "Pi(raw)", "Pi(K80)", "Pi(Trs)", "Pi(Trv)", "Pi(Syn)", "Pi(Nsyn)", "Pi(dN2ds)",
                      "Ratio(raw)", "Ratio(K80)", "Ratio(Trs)", "Ratio(Trv)", "Ratio(Syn)", "Ratio(Nsyn)", "Ratio(dN2ds)",
                      "FAD(raw)", "FAD(K80)", "FAD(Trs)", "FAD(Trv)", "FAD(Syn)", "FAD(Nsyn)", "FAD(dN2ds)")
for(i in 1:length(subjects.placebo)){
  sub <- subjects.placebo[i]
  data <- TEP.placebo %>% filter(subject==sub) %>% select(nuc_sequence, reads)
  seq <- as.vector(data$nuc_sequence)
  count <- as.vector(data$reads)
  nm <- SortByMutations(seq,count)$nm
  len <- nchar(seq[1])
  dst <- DNA.dist(seq)
  placebo.df[i,1:43] <- c(Shannon.entropy(count), N.Shannon.entropy(count), 
                      GiniSimpson(count), 
                      qD(count,0), qD(count,1), qD(count,2), qD(count,3), 
                      renyi(count,0), renyi(count,2), HCq(count,3),
                      poly.sites(seq)$psites, poly.sites(seq)$nmuts, 
                      MutationFreq(SortByMutations(seq,count)$nr,SortByMutations(seq,count)$nm,len), correct.seg.sites(seq, count),
                      rao(count,dst), 
                      sq.diversity2(seq,count,1)$Mf[1], sq.diversity2(seq,count,1)$Mf[2], sq.diversity2(seq,count,1)$Mf[3], sq.diversity2(seq,count,1)$Mf[4], 
                      sq.diversity2(seq,count,1)$Mf[5], sq.diversity2(seq,count,1)$Mf[6], sq.diversity2(seq,count,1)$Mf[7],
                      sq.diversity2(seq,count,1)$Pi[1], sq.diversity2(seq,count,1)$Pi[2], sq.diversity2(seq,count,1)$Pi[3], sq.diversity2(seq,count,1)$Pi[4],
                      sq.diversity2(seq,count,1)$Pi[5], sq.diversity2(seq,count,1)$Pi[6], sq.diversity2(seq,count,1)$Pi[7],
                      sq.diversity2(seq,count,1)$Ratio[1], sq.diversity2(seq,count,1)$Ratio[2], sq.diversity2(seq,count,1)$Ratio[3], sq.diversity2(seq,count,1)$Ratio[4],
                      sq.diversity2(seq,count,1)$Ratio[5], sq.diversity2(seq,count,1)$Ratio[6], sq.diversity2(seq,count,1)$Ratio[7],
                      sq.diversity2(seq,count,1)$FAD[1], sq.diversity2(seq,count,1)$FAD[2], sq.diversity2(seq,count,1)$FAD[3], sq.diversity2(seq,count,1)$FAD[4],
                      sq.diversity2(seq,count,1)$FAD[5], sq.diversity2(seq,count,1)$FAD[6], sq.diversity2(seq,count,1)$FAD[7])
}
rownames(placebo.df) <- subjects.placebo

#save as csv file
write.csv(placebo.df, "C:/Users/Swu2/Documents/placebo.df")


symbox(~Mfmax, data=placebo.df, start=0.0001)
aa.nms <- c("A","R","N","D", "B", "C", "E","Q","Z",
            "G","H","I","L","K","M","F","P","S","T",
            "W","Y","V")

expland.grid(g1,g2,g3,g4,g5,g6,g7,g8)


#=========================================AA analysis
source("C:/Users/Swu2/Documents/AADiversityIndices.R")
setwd("T:/vaccine/rtss_malaria_sieve/")
TEP.data <- read_tsv("qdata/sequences/TEP.tsv")

#subset for subjects age 5-17 months, separated into placebo and vaccine
clinical.data <- read_csv("adata/RTSSclinicalData.csv")
clinical.data <- subset(clinical.data, ageCateg=="[5-17] months")
clinical.placebo <- subset(clinical.data, vaccine==0)
clinical.vaccine <- subset(clinical.data, vaccine==1)

#TEP data for all participants age 5-17 months, as well as separated into vaccine and placebo groups
TEP.all <- TEP.data[TEP.data$subject %in% clinical.data$id,]
TEP.placebo <- TEP.data[TEP.data$subject %in% clinical.placebo$id,]
TEP.vaccine <- TEP.data[TEP.data$subject %in% clinical.vaccine$id,]

#number of distinct subjects in TEP data for participants age 5-17 months
subjects.TEP <- TEP.all %>% distinct(subject)
subjects.TEP <- subjects.TEP$subject

#number of distinct subjects in TEP placebo group for participants age 5-17 months
subjects.placebo <- TEP.placebo %>% distinct(subject)
subjects.placebo <- subjects.placebo$subject

#=======================TEP AA dataframe of indices================================

TEP.aa.df <- data.frame(matrix(ncol = 17, nrow = length(subjects.TEP)))
colnames(TEP.aa.df) <- c("Sh", "N.Sh", 
                      "GS",
                      "H(0)","H(1)","H(2)", "H(3)",
                      "Ren(0)","Ren(2)","HC(3)",
                      "PolySites","UniqMuts", "CorSegSites",
                      "Mf",  
                      "Rao",
                      "Ratio",
                      "Hamming")
for(i in 1:length(subjects.TEP)){
  sub <- subjects.TEP[i]
  data <- TEP.all %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  len <- nchar(seq[1])
  dst <- stringDist(seq)/nchar(seq[1])
  lst <- Aa.poly.sites(seq)
  TEP.aa.df[i,1:17] <- c(Shannon.entropy(count), N.Shannon.entropy(count),
                      GiniSimpson(count),
                      qD(count,0), qD(count,1), qD(count,2), qD(count,3),
                      renyi(count,0), renyi(count,2), HCq(count,3),
                      Aa.poly.sites(seq)$psites, Aa.poly.sites(seq)$nmuts, AA.correct.seg.sites(seq, count),
                      Aa.sq.diversity(seq,count,nm)$Mf, Aa.sq.diversity(seq,count,nm)$Rao,
                      Aa.sq.diversity(seq,count,nm)$Ratio, 
                      mean(nm))
}
rownames(TEP.aa.df) <- subjects.TEP
vaccine <- ifelse(rownames(TEP.aa.df) %in% TEP.vaccine$subject, "vaccine", "placebo")
#save as csv file
write.csv(TEP.aa.df, "C:/Users/Swu2/Documents/TEP.aa.df")


#=========================Placebo AA
placebo.aa.df <- data.frame(matrix(ncol = 17, nrow = length(subjects.placebo)))
colnames(placebo.aa.df) <- c("Sh", "N.Sh", 
                         "GS",
                         "H(0)","H(1)","H(2)", "H(3)",
                         "Ren(0)","Ren(2)","HC(3)",
                         "PolySites","UniqMuts", "CorSegSites",
                         "Mf",  
                         "Rao",
                         "Ratio",
                         "Hamming")
for(i in 1:length(subjects.placebo)){
  sub <- subjects.placebo[i]
  data <- TEP.placebo %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  len <- nchar(seq[1])
  dst <- stringDist(seq)/nchar(seq[1])
  lst <- Aa.poly.sites(seq)
  placebo.aa.df[i,1:17] <- c(Shannon.entropy(count), N.Shannon.entropy(count), 
                         GiniSimpson(count), 
                         qD(count,0), qD(count,1), qD(count,2), qD(count,3), 
                         renyi(count,0), renyi(count,2), HCq(count,3),
                         Aa.poly.sites(seq)$psites, Aa.poly.sites(seq)$nmuts, AA.correct.seg.sites(seq, count),
                         Aa.sq.diversity(seq,count,nm)$Mf, Aa.sq.diversity(seq,count,nm)$Rao, 
                         Aa.sq.diversity(seq,count,nm)$Ratio,
                         mean(nm))
}
rownames(placebo.aa.df) <- subjects.placebo

#impute infinite values with 0
is.na(placebo.aa.df) <- do.call(cbind,lapply(placebo.aa.df, is.infinite))
placebo.aa.df[is.na(placebo.aa.df)] <- 0

#save as csv file
write.csv(placebo.aa.df, "C:/Users/Swu2/Documents/placebo.aa.df")



#============================create fringe-trimmed aa placebo dataframe
placebo.aa.fringe <- data.frame(matrix(ncol = 17, nrow = length(subjects.placebo)))
colnames(placebo.aa.fringe) <- c("Sh", "N.Sh", 
                                 "GS",
                                 "H(0)","H(1)","H(2)", "H(3)",
                                 "Ren(0)","Ren(2)","HC(3)",
                                 "PolySites","UniqMuts", "CorSegSites",
                                 "Mf",  
                                 "Rao",
                                 "Ratio",
                                 "Hamming")

for(i in 1:length(subjects.placebo)){
  sub <- subjects.placebo[i]
  data <- TEP.placebo %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  len <- nchar(seq[1])
  dst <- stringDist(seq)/nchar(seq[1])
  lst <- Aa.poly.sites(seq)
  placebo.aa.fringe[i,1:17] <- c(Shannon.entropy(count), N.Shannon.entropy(count), 
                                 GiniSimpson(count), 
                                 qD(count,0), qD(count,1), qD(count,2), qD(count,3), 
                                 renyi(count,0), renyi(count,2), HCq(count,3),
                                 Aa.poly.sites(seq)$psites, Aa.poly.sites(seq)$nmuts, AA.correct.seg.sites(seq, count),
                                 Aa.sq.diversity(seq,count,nm)$Mf, Aa.sq.diversity(seq,count,nm)$Rao, 
                                 Aa.sq.diversity(seq,count,nm)$Ratio,
                                 mean(nm))
}
rownames(placebo.aa.fringe) <- subjects.placebo

#impute infinite values with 0
is.na(placebo.aa.fringe) <- do.call(cbind,lapply(placebo.aa.fringe, is.infinite))
placebo.aa.fringe[is.na(placebo.aa.fringe)] <- 0

#save as csv file
write.csv(placebo.aa.fringe, "C:/Users/Swu2/Documents/placebo.aa.fringe")



#=================================create fringe-trimmed aa vaccine dataframe
vaccine.aa.fringe <- data.frame(matrix(ncol = 17, nrow = length(subjects.vaccine)))
colnames(vaccine.aa.fringe) <- c("Sh", "N.Sh", 
                                 "GS",
                                 "H(0)","H(1)","H(2)", "H(3)",
                                 "Ren(0)","Ren(2)","HC(3)",
                                 "PolySites","UniqMuts", "CorSegSites",
                                 "Mf",  
                                 "Rao",
                                 "Ratio",
                                 "Hamming")

for(i in 1:length(subjects.vaccine)){
  sub <- subjects.vaccine[i]
  data <- TEP.vaccine %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  len <- nchar(seq[1])
  dst <- stringDist(seq)/nchar(seq[1])
  lst <- Aa.poly.sites(seq)
  vaccine.aa.fringe[i,1:17] <- c(Shannon.entropy(count), N.Shannon.entropy(count), 
                                 GiniSimpson(count), 
                                 qD(count,0), qD(count,1), qD(count,2), qD(count,3), 
                                 renyi(count,0), renyi(count,2), HCq(count,3),
                                 Aa.poly.sites(seq)$psites, Aa.poly.sites(seq)$nmuts, AA.correct.seg.sites(seq, count),
                                 Aa.sq.diversity(seq,count,nm)$Mf, Aa.sq.diversity(seq,count,nm)$Rao, 
                                 Aa.sq.diversity(seq,count,nm)$Ratio,
                                 mean(nm))
}
rownames(vaccine.aa.fringe) <- subjects.vaccine

#impute infinite values with 0
is.na(vaccine.aa.fringe) <- do.call(cbind,lapply(vaccine.aa.fringe, is.infinite))
vaccine.aa.fringe[is.na(vaccine.aa.fringe)] <- 0

#save as csv file
write.csv(vaccine.aa.fringe, "C:/Users/Swu2/Documents/vaccine.aa.fringe")

#==================================Clinical/parasite comparisons==========================
#create dataframe of clinical vs. parasite pos comparisons for select indices
comparisons.df <- data.frame(matrix(ncol = 16, nrow = length(subjects.both)))
colnames(comparisons.df) <- c("H0.c", "H1.c","H2.c","CorSegSites.c","GS.c","Pi.c","Mf.c","Hamming.c",
                              "H0.x", "H1.x","H2.x","CorSegSites.x","GS.x","Pi.x","Mf.x","Hamming.x")

for(i in 1:length(subjects.both)){
  sub <- subjects.both[i]
  #clinical cases
  data <- clinical %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comparisons.df[i,1:8] <- c(qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                             GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                             mean(nm))
  #parasite positive cases
  data <- parasitepos %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comparisons.df[i,9:16] <- c(qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                              GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                              mean(nm))
}
rownames(comparisons.df) <- subjects.both

#impute infinite values with 0
is.na(comparisons.df) <- do.call(cbind,lapply(comparisons.df, is.infinite))
comparisons.df[is.na(comparisons.df)] <- 0

#save as csv file
write.csv(comparisons.df, "C:/Users/Swu2/Documents/ClinicalParasiteComp")

#========================Clnical/Parasite for all AA Indices=====================

comps <- data.frame(matrix(ncol = 30, nrow = length(subjects.both)))
colnames(comps) <- c("Sh.c", "N.Sh.c", "Ren0.c","Ren2.c","HC3.c",
                     "PolySites.c","UniqMuts.c","H0.c", "H1.c","H2.c","CorSegSites.c","GS.c",
                     "Pi.c","Mf.c","Hamming.c",
                     "Sh.x", "N.Sh.x", "Ren0.x","Ren2.x","HC3.x",
                     "PolySites.x","UniqMuts.x", "H0.x", "H1.x","H2.x","CorSegSites.x","GS.x",
                     "Pi.x","Mf.x","Hamming.x")
for(i in 1:length(subjects.both)){
  sub <- subjects.both[i]
  #clinical cases
  data <- clinical %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comps[i,1:15] <- c(Shannon.entropy(count), N.Shannon.entropy(count),  
                     renyi(count,0), renyi(count,2), HCq(count,3),
                     Aa.poly.sites(seq)$psites, Aa.poly.sites(seq)$nmuts, 
                     qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                     GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                     mean(nm))
  #parasite positive cases
  data <- parasitepos %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comps[i,16:30] <- c(Shannon.entropy(count), N.Shannon.entropy(count),  
                      renyi(count,0), renyi(count,2), HCq(count,3),
                      Aa.poly.sites(seq)$psites, Aa.poly.sites(seq)$nmuts, 
                      qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                      GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                      mean(nm))
}
rownames(comps) <- subjects.both

#impute infinite values with 0
is.na(comps) <- do.call(cbind,lapply(comps, is.infinite))
comps[is.na(comps)] <- 0

#save as csv file
write.csv(comps, "C:/Users/Swu2/Documents/Comps")

#==================================Clinical/parasite comparisons for vaccinees==========================
#create dataframe of clinical vs. parasite pos comparisons for select indices
comparisons.df <- data.frame(matrix(ncol = 16, nrow = length(subjects.both)))
colnames(comparisons.df) <- c("H0.c", "H1.c","H2.c","CorSegSites.c","GS.c","Pi.c","Mf.c","Hamming.c",
                              "H0.x", "H1.x","H2.x","CorSegSites.x","GS.x","Pi.x","Mf.x","Hamming.x")

for(i in 1:length(subjects.both)){
  sub <- subjects.both[i]
  #clinical cases
  data <- clinical %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comparisons.df[i,1:8] <- c(qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                             GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                             mean(nm))
  #parasite positive cases
  data <- parasitepos %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comparisons.df[i,9:16] <- c(qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                              GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                              mean(nm))
}
rownames(comparisons.df) <- subjects.both

#impute infinite values with 0
is.na(comparisons.df) <- do.call(cbind,lapply(comparisons.df, is.infinite))
comparisons.df[is.na(comparisons.df)] <- 0

#save as csv file
write.csv(comparisons.df, "C:/Users/Swu2/Documents/ClinicalParasiteCompVax")


#==================================Clinical/parasite comparisons for SERA2 Placebos==========================
#create dataframe of clinical vs. parasite pos comparisons for select indices
comparisons.df <- data.frame(matrix(ncol = 16, nrow = length(subjects.both)))
colnames(comparisons.df) <- c("H0.c", "H1.c","H2.c","CorSegSites.c","GS.c","Pi.c","Mf.c","Hamming.c",
                              "H0.x", "H1.x","H2.x","CorSegSites.x","GS.x","Pi.x","Mf.x","Hamming.x")

for(i in 1:length(subjects.both)){
  sub <- subjects.both[i]
  #clinical cases
  data <- clinical %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comparisons.df[i,1:8] <- c(qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                             GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                             mean(nm))
  #parasite positive cases
  data <- parasitepos %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comparisons.df[i,9:16] <- c(qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                              GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                              mean(nm))
}
rownames(comparisons.df) <- subjects.both

#impute infinite values with 0
is.na(comparisons.df) <- do.call(cbind,lapply(comparisons.df, is.infinite))
comparisons.df[is.na(comparisons.df)] <- 0

#save as csv file
write.csv(comparisons.df, "C:/Users/Swu2/Documents/ClinicalParasiteCompSERA2")


#==================================Clinical/parasite comparisons for SERA2 Vaccinees==========================
#create dataframe of clinical vs. parasite pos comparisons for select indices
comparisons.df <- data.frame(matrix(ncol = 16, nrow = length(subjects.both)))
colnames(comparisons.df) <- c("H0.c", "H1.c","H2.c","CorSegSites.c","GS.c","Pi.c","Mf.c","Hamming.c",
                              "H0.x", "H1.x","H2.x","CorSegSites.x","GS.x","Pi.x","Mf.x","Hamming.x")

for(i in 1:length(subjects.both)){
  sub <- subjects.both[i]
  #clinical cases
  data <- clinical %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comparisons.df[i,1:8] <- c(qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                             GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                             mean(nm))
  #parasite positive cases
  data <- parasitepos %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  #fringe-trimming
  cutoff <- qbinom(0.9, size=sum(count), prob = 0.01, lower.tail=TRUE)
  for(j in 1:length(count)){
    if(count[j] < cutoff){
      seq[j] <- NA
      count[j] <- NA
    }
  }
  seq <- na.omit(seq)
  count <- na.omit(count)
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=vax3D7)
  nm <- nmismatch(psa)
  #calculate indices
  comparisons.df[i,9:16] <- c(qD(count,0), qD(count,1), qD(count,2), AA.correct.seg.sites(seq, count),
                              GiniSimpson(count), Aa.sq.diversity(seq,count,nm)$Rao, Aa.sq.diversity(seq,count,nm)$Mf, 
                              mean(nm))
}
rownames(comparisons.df) <- subjects.both

#impute infinite values with 0
is.na(comparisons.df) <- do.call(cbind,lapply(comparisons.df, is.infinite))
comparisons.df[is.na(comparisons.df)] <- 0

#save as csv file
write.csv(comparisons.df, "C:/Users/Swu2/Documents/ClinicalParasiteCompSERA2Vax")

#===========================================================



which(subjects.TEP==2177)
for(i in 1:length(subjects.TEP)){
  sub <- subjects.TEP[1470]
  data <- TEP.all %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  master <- seq[which.max(count)]
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=master)
  nm <- nmismatch(psa)
  len <- nchar(seq[1])
  dst <- stringDist(seq)/nchar(seq[1])
  lst <- Aa.poly.sites(seq)
  TEP.test[i] <- Aa.sq.diversity(seq, count, nm)
}

for(i in 1:length(subjects.TEP)){
  sub <- subjects.TEP[3]
  data <- TEP.all %>% filter(subject==sub) %>% select(pep_sequence, reads)
  seq <- as.vector(data$pep_sequence)
  count <- as.vector(data$reads)
  dst <- stringDist(seq)/nchar(seq[1])
  rao.aa(count, dst)
  
  master <- seq[which.max(count)]
  bseqs <- AAStringSet(seq)
  psa <- pairwiseAlignment(pattern=bseqs,subject=master)
  nm <- nmismatch(psa)
  MutationFreq(count,nm,nchar(seq[1]))
}

#haplotype and reads summary
num.haplotypes <- numeric(length(subjects.placebo))
num.reads <- numeric(length(subjects.placebo))
for(i in 1:length(subjects.placebo)){
  sub <- subjects.TEP[i]
  data <- TEP.all %>% filter(subject==sub) %>% select(nuc_sequence, reads)
  count <- as.vector(data$reads)
  num.haplotypes[i] <- length(count)
  num.reads[i] <- sum(count)
}
summary(num.haplotypes)
summary(num.reads)
count.avg <- mean(num.haplotypes)
reads.avg <- mean(num.reads)