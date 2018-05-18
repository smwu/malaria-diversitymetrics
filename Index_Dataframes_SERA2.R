#vaccine 3D7 strain
vax3D7 <- "LSSDGSRVTTQARIEKPKQQPTLPTLAQETQPQQQQQQKEVGSGIGAEQKVESARPGAEVSQSDVERAGRSSGTGGSVGTKISPG"  #85 bases

#load SERA2 data and clinical and parasite indicator files
setwd("T:/vaccine/rtss_malaria_sieve/")
SERA2.data <- read_tsv("qdata/sequences/SERA2.tsv")
marks.data.x <- read_tsv("adata/marks_data_x_sites.tsv")  #parasite positive cases
marks.data.c <- read_tsv("adata/marks_data_c_sites.tsv")  #clinical malaria cases

#subset for subjects age 5-17 months, separated into placebo and vaccine
clinical.data <- read_csv("adata/RTSSclinicalData.csv")
clinical.data <- subset(clinical.data, ageCateg=="[5-17] months")
clinical.placebo <- subset(clinical.data, vaccine==0)
clinical.vaccine <- subset(clinical.data, vaccine==1)

#SERA2 data for all participants age 5-17 months, as well as separated into vaccine and placebo groups
SERA2.all <- SERA2.data[SERA2.data$subject %in% clinical.data$id,]
SERA2.placebo <- SERA2.data[SERA2.data$subject %in% clinical.placebo$id,]
SERA2.vaccine <- SERA2.data[SERA2.data$subject %in% clinical.vaccine$id,]

#number of distinct subjects in SERA2 data for participants age 5-17 months
subjects.SERA2 <- SERA2.all %>% distinct(subject)
subjects.SERA2 <- subjects.SERA2$subject

#number of distinct subjects in SERA2 placebo group for participants age 5-17 months
subjects.placebo <- SERA2.placebo %>% distinct(subject)
subjects.placebo <- subjects.placebo$subject

#number of distinct subjects in SERA2 vaccine group for participants age 5-17 months
subjects.vaccine <- SERA2.vaccine %>% distinct(subject)
subjects.vaccine <- subjects.vaccine$subject



#==================================Clinical/parasite comparisons for SERA2 Placebos==========================

#cases of clinical malaria for SERA2 placebo subjects aged 5-17 months
clinical <- SERA2.data[SERA2.data$subject %in% subjects.placebo,]
clinical <- clinical[clinical$sample %in% marks.data.c$sample,]

#cases of parasite positivity for SERA2 placebo subjects aged 5-17 months
parasitepos <- SERA2.data[SERA2.data$subject %in% subjects.placebo,]
parasitepos <- parasitepos[parasitepos$sample %in% marks.data.x$sample,]

#cases of both clinical malaria and parasite positivity for SERA2 placebo subjects aged 5-17 months
both <- clinical[clinical$subject %in% parasitepos$subject,]
both <- both[,1:11]

#number of distinct subjects with both clinical and parasite endpoints 
#for SERA2 placebo subjects aged 5-17 months
subjects.both <- both %>% distinct(subject)
subjects.both <- subjects.both$subject


#create dataframe of clinical vs. parasite pos comparisons for select indices
comparisons.df <- data.frame(matrix(ncol = 16, nrow = length(subjects.both)))
colnames(comparisons.df) <- c("H0.c", "H1.c","H2.c","CorSegSites.c","PiWt.c","PiEq.c","HamWt.c","HamEq.c",
                              "H0.x", "H1.x","H2.x","CorSegSites.x","PiWt.x","PiEq.x","HamWt.x","HamEq.x")

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
  dst <- stringDist(seq)
  #calculate indices
  comparisons.df[i,1:8] <- c(qD(count,0), qD(count,1), qD(count,2), 
                              AA.correct.seg.sites(seq, count),
                              rao.aa(count, dst), pi.aa(count,dst), 
                              sum(nm*count)/sum(count), mean(nm))
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
  dst <- stringDist(seq)
  #calculate indices
  comparisons.df[i,9:16] <- c(qD(count,0), qD(count,1), qD(count,2), 
                              AA.correct.seg.sites(seq, count),
                              rao.aa(count, dst), pi.aa(count,dst), 
                              sum(nm*count)/sum(count), mean(nm))
}
rownames(comparisons.df) <- subjects.both

#impute infinite values with 0
is.na(comparisons.df) <- do.call(cbind,lapply(comparisons.df, is.infinite))
comparisons.df[is.na(comparisons.df)] <- 0

#save as csv file
write.csv(comparisons.df, "C:/Users/Swu2/Documents/ClinicalParasiteCompSERA2")



#==================================Clinical/parasite comparisons for SERA2 Vaccinees==========================

#cases of clinical malaria for SERA2 placebo subjects aged 5-17 months
clinical <- SERA2.data[SERA2.data$subject %in% subjects.vaccine,]
clinical <- clinical[clinical$sample %in% marks.data.c$sample,]

#cases of parasite positivity for SERA2 placebo subjects aged 5-17 months
parasitepos <- SERA2.data[SERA2.data$subject %in% subjects.vaccine,]
parasitepos <- parasitepos[parasitepos$sample %in% marks.data.x$sample,]

#cases of both clinical malaria and parasite positivity for SERA2 placebo subjects aged 5-17 months
both <- clinical[clinical$subject %in% parasitepos$subject,]
both <- both[,1:11]

#number of distinct subjects with both clinical and parasite endpoints 
#for SERA2 placebo subjects aged 5-17 months
subjects.both <- both %>% distinct(subject)
subjects.both <- subjects.both$subject


#create dataframe of clinical vs. parasite pos comparisons for select indices
comparisons.df <- data.frame(matrix(ncol = 16, nrow = length(subjects.both)))
colnames(comparisons.df) <- c("H0.c", "H1.c","H2.c","CorSegSites.c","PiWt.c","PiEq.c","HamWt.c","HamEq.c",
                              "H0.x", "H1.x","H2.x","CorSegSites.x","PiWt.x","PiEq.x","HamWt.x","HamEq.x")

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
  dst <- stringDist(seq)
  #calculate indices
  comparisons.df[i,1:8] <- c(qD(count,0), qD(count,1), qD(count,2), 
                              AA.correct.seg.sites(seq, count),
                              rao.aa(count, dst), pi.aa(count,dst), 
                              sum(nm*count)/sum(count), mean(nm))
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
  dst <- stringDist(seq)
  #calculate indices
  comparisons.df[i,9:16] <- c(qD(count,0), qD(count,1), qD(count,2), 
                              AA.correct.seg.sites(seq, count),
                              rao.aa(count, dst), pi.aa(count,dst), 
                              sum(nm*count)/sum(count), mean(nm))
}
rownames(comparisons.df) <- subjects.both

#impute infinite values with 0
is.na(comparisons.df) <- do.call(cbind,lapply(comparisons.df, is.infinite))
comparisons.df[is.na(comparisons.df)] <- 0

#save as csv file
write.csv(comparisons.df, "C:/Users/Swu2/Documents/ClinicalParasiteCompSERA2Vax")