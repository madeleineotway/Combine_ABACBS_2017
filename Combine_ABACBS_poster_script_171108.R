######################################
##         Madeleine J Otway        ##
##        motway@cmri.org.au        ##
##            08/11/2017            ##
##     Combine and ABACBS 2017      ##
##      Source code for poster      ##
##     Data available on GitHub     ##
##   R version 3.4.0 (2017-04-21)   ##
##     Bioconductor version 3.6     ##
######################################

#install.source()
library(dialects)

#install.packages("fBasics")
library(fBasics)

#install.packages("ggplot2")
library(ggplot2)

#source("http://bioconductor.org/biocLite.R")
#biocLite("plyr")
library(plyr)

##
## Human SRL convert to Rat
human.srl <- import.srl("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Human_SRL_dialects.txt")
rat.fasta <- import.fasta("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/Fasta/uniprot-rat.fasta")
digest.rat <- digest.fasta(rat.fasta)
rat.from.human <- convert.species(human.srl, digest.rat)
export.srl(rat.from.human, "/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Rat_from_Human_SRL_dialects.txt")

##

## Rat SRL convert to Human
rat.srl <- import.srl("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Rat_SRL_dialects.txt")
human.fasta <- import.fasta("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/Fasta/uniprot-human.fasta")
digest.human <- digest.fasta(human.fasta)
human.from.rat <- convert.species(rat.srl, digest.human)
export.srl(human.from.rat, "/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Human_from_Rat_SRL_dialects.txt")


##


## Import results from rat SRL
pep.rat <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Rat/MJO_PGH_Rat_PeakExtraction_dialects_areas - Peptides.txt",
                      sep = "\t", header = T)
fdr.rat <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Rat/MJO_PGH_Rat_PeakExtraction_dialects_FDR.txt",
                      sep = "\t", header = T)
## Remove decoy and match FDR and peptide data frames
fdr.decoy <- fdr.rat[!(fdr.rat$Decoy=="True"),]
fdr.decoy$fdr.pep <- ifelse(fdr.decoy$Protein %in% pep.rat$Protein &
                              fdr.decoy$Peptide %in% pep.rat$Peptide &
                              fdr.decoy$Precursor.MZ %in% pep.rat$Precursor.MZ,
                            1, 0)
fdr.match <- fdr.decoy[!(fdr.decoy$fdr.pep == "0"),]
## Remove all peptides with an FDR greater that 1%
sort.fdr <- arrange(fdr.match, Protein, Peptide, Precursor.MZ)
sort.pep <- arrange(pep.rat, Protein, Peptide, Precursor.MZ)
sort.fdr$Label <- NULL
sort.fdr$Decoy <- NULL
rat.01 <- sort.pep
for (i in 6:10) {
  rat.01[,i] <- ifelse(sort.fdr[,i] <= 0.01, rat.01[,i], NA)
}
## Remove Missing values
rat.01 <- rat.01[(rowSums(is.na(rat.01[6:10])) <= 0) ,]

##

## Repeat for Rat SRL from Human SRL results
pep.rat.h <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Rat/MJO_PGH_Rat_from_Human_PeakExtraction_dialects_areas - Peptides.txt",
                        sep = "\t", header = T)
fdr.rat.h <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Rat/MJO_PGH_Rat_from_Human_PeakExtraction_dialects_FDR.txt",
                        sep = "\t", header = T)
fdr.decoy.h <- fdr.rat.h[!(fdr.rat.h$Decoy=="True"),]
fdr.decoy.h$fdr.pep <- ifelse(fdr.decoy.h$Protein %in% pep.rat.h$Protein &
                                fdr.decoy.h$Peptide %in% pep.rat.h$Peptide &
                                fdr.decoy.h$Precursor.MZ %in% pep.rat.h$Precursor.MZ,
                              1, 0)
fdr.match.h <- fdr.decoy.h[!(fdr.decoy.h$fdr.pep == "0"),]
sort.fdr.h <- arrange(fdr.match.h, Protein, Peptide, Precursor.MZ)
sort.pep.h <- arrange(pep.rat.h, Protein, Peptide, Precursor.MZ)
sort.fdr.h$Label <- NULL
sort.fdr.h$Decoy <- NULL
rat.h.01 <- sort.pep.h
for (i in 6:10) {
  rat.h.01[,i] <- ifelse(sort.fdr.h[,i] <= 0.01, rat.h.01[,i], NA)
}
rat.h.01 <- rat.h.01[(rowSums(is.na(rat.h.01[6:10])) <= 0) ,]

##

## Import results from rat SRL
pep.human <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Human/MJO_PGH_Human_PeakExtraction_dialects_areas - Peptides.txt",
                        sep = "\t", header = T)
fdr.human <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Human/MJO_PGH_Human_PeakExtraction_dialects_FDR.txt",
                        sep = "\t", header = T)
## Remove decoy and match FDR and peptide data frames
fdr.decoy <- fdr.human[!(fdr.human$Decoy=="True"),]
fdr.decoy$fdr.pep <- ifelse(fdr.decoy$Protein %in% pep.human$Protein &
                              fdr.decoy$Peptide %in% pep.human$Peptide &
                              fdr.decoy$Precursor.MZ %in% pep.human$Precursor.MZ,
                            1, 0)
fdr.match <- fdr.decoy[!(fdr.decoy$fdr.pep == "0"),]
## Remove all peptides with an FDR greater that 1%
sort.fdr <- arrange(fdr.match, Protein, Peptide, Precursor.MZ)
sort.pep <- arrange(pep.human, Protein, Peptide, Precursor.MZ)
sort.fdr$Label <- NULL
sort.fdr$Decoy <- NULL
human.01 <- sort.pep
for (i in 6:10) {
  human.01[,i] <- ifelse(sort.fdr[,i] <= 0.01, human.01[,i], NA)
}
## Remove Missing values
human.01 <- human.01[(rowSums(is.na(human.01[6:10])) <= 0) ,]

##
## Repeat for Rat SRL from Human SRL results
pep.human.r <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Human/MJO_PGH_Human_from_Rat_PeakExtraction_dialects_areas - Peptides.txt",
                          sep = "\t", header = T)
fdr.human.r <- read.table("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/PeakExtractions/Human/MJO_PGH_Human_from_Rat_PeakExtraction_dialects_FDR.txt",
                          sep = "\t", header = T)
fdr.decoy.h <- fdr.human.r[!(fdr.human.r$Decoy=="True"),]
fdr.decoy.h$fdr.pep <- ifelse(fdr.decoy.h$Protein %in% pep.human.r$Protein &
                                fdr.decoy.h$Peptide %in% pep.human.r$Peptide &
                                fdr.decoy.h$Precursor.MZ %in% pep.human.r$Precursor.MZ,
                              1, 0)
fdr.match.h <- fdr.decoy.h[!(fdr.decoy.h$fdr.pep == "0"),]
sort.fdr.h <- arrange(fdr.match.h, Protein, Peptide, Precursor.MZ)
sort.pep.h <- arrange(pep.human.r, Protein, Peptide, Precursor.MZ)
sort.fdr.h$Label <- NULL
sort.fdr.h$Decoy <- NULL
human.r.01 <- sort.pep.h
for (i in 6:10) {
  human.r.01[,i] <- ifelse(sort.fdr.h[,i] <= 0.01, human.r.01[,i], NA)
}
human.r.01 <- human.r.01[(rowSums(is.na(human.r.01[6:10])) <= 0) ,]

##
#install.packages("fBasics")
#library(fBasics)


human.01$human.match <- ifelse(human.01$Protein %in% human.r.01$Protein &
                                 human.01$Peptide %in% human.r.01$Peptide &
                                 human.01$Precursor.MZ %in% human.r.01$Precursor.MZ &
                                 human.01$Precursor.Charge %in% human.r.01$Precursor.Charge,
                               1, 0)
human.match <- human.01[!(human.01$human.match == "0"),]

human.r.01$human.match <- ifelse(human.r.01$Protein %in% human.01$Protein &
                                   human.r.01$Peptide %in% human.01$Peptide &
                                   human.r.01$Precursor.MZ %in% human.01$Precursor.MZ &
                                   human.r.01$Precursor.Charge %in% human.01$Precursor.Charge,
                                 1, 0)

human.r.match <- human.r.01[!(human.r.01$human.match == "0"),]


human.r.match$mean <- rowMeans(human.r.match[4:10])
human.r.match$sd <- rowStdevs(human.r.match[4:10])
human.match$mean <- rowMeans(human.match[4:10])
human.match$sd <- rowStdevs(human.match[4:10])

rat.srl <- import.srl("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Rat_SRL_dialects.txt")
rat.srl$rat.match <- ifelse(rat.srl$stripped_sequence %in% human.match$Peptide &
                              rat.srl$prec_z %in% human.match$Precursor.Charge,
                            1, 0)
rat.srl.match <- rat.srl[!(rat.srl$rat.match == "0"),]
rat.srl.match$new.u <- rat.match$uniprot_id
rat.srl.match$new.u <- gsub("..\\|.*\\|", "", rat.match$new.u)
rat.srl.match$new.u <- gsub("_RAT", "", rat.match$new.u)

human.match$new.u <- human.match$Protein
human.match$new.u <- gsub("..\\|.*\\|", "", human.match$new.u)
human.match$new.u <- gsub("_HUMAN", "", human.match$new.u)
human.match$u.match <- ifelse(human.match$new.u %in% rat.srl.match$new.u, 
                              1, 0)
human.u.match <- human.match[!(human.match$u.match == "1"),]
human.u.match <- human.u.match[!duplicated(human.u.match$new.u),]


human.r.match$new.u <- human.r.match$Protein
human.r.match$new.u <- gsub("..\\|.*\\|", "", human.r.match$new.u)
human.r.match$new.u <- gsub("_HUMAN", "", human.r.match$new.u)
human.r.match$u.match <- ifelse(human.r.match$new.u %in% rat.srl.match$new.u, 
                                1, 0)
human.r.u.match <- human.r.match[!(human.r.match$u.match == "1"),]
human.r.u.match <- human.r.u.match[!duplicated(human.r.u.match$new.u),]


h <- data.frame("Peptide" = human.u.match[1:5,2])
h$Intensity <- human.u.match[1:5,12]
h$group <- "human"

hr <- data.frame("Peptide" = human.r.u.match[1:5,2])
hr$Intensity <- human.r.u.match[1:5,12]
hr$group <- "human.r"
h.hr <- rbind(h, hr)

ggplot(h.hr, aes(x = Peptide, y = Intensity, fill = group)) +
  geom_bar(stat = "identity", color = "black", position=position_dodge()) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 20),
        axis.text.y = element_text(hjust = 1, size = 20)) +
  guides(fill = F) +
  scale_fill_manual(values=c('#008000','#FF0080'))

##

##
rat.01$rat.match <- ifelse(rat.01$Protein %in% rat.h.01$Protein &
                             rat.01$Peptide %in% rat.h.01$Peptide &
                             rat.01$Precursor.MZ %in% rat.h.01$Precursor.MZ &
                             rat.01$Precursor.Charge %in% rat.h.01$Precursor.Charge,
                           1, 0)
rat.match <- rat.01[!(rat.01$rat.match == "0"),]

rat.h.01$rat.match <- ifelse(rat.h.01$Protein %in% rat.01$Protein &
                               rat.h.01$Peptide %in% rat.01$Peptide &
                               rat.h.01$Precursor.MZ %in% rat.01$Precursor.MZ &
                               rat.h.01$Precursor.Charge %in% rat.01$Precursor.Charge,
                             1, 0)

rat.h.match <- rat.h.01[!(rat.h.01$rat.match == "0"),]


rat.h.match$mean <- rowMeans(rat.h.match[4:10])
rat.h.match$sd <- rowStdevs(rat.h.match[4:10])
rat.match$mean <- rowMeans(rat.match[4:10])
rat.match$sd <- rowStdevs(rat.match[4:10])

human.srl <- import.srl("/Users/MadeleineOtway/Documents/PhD/R/dialects/Data/SRL/MJO_PGH_Human_SRL_dialects.txt")
human.srl$rat.match <- ifelse(human.srl$stripped_sequence %in% rat.match$Peptide &
                                human.srl$prec_z %in% rat.match$Precursor.Charge,
                              1, 0)
human.match <- human.srl[!(human.srl$rat.match == "0"),]
human.match$new.u <- human.match$uniprot_id
human.match$new.u <- gsub("..\\|.*\\|", "", human.match$new.u)
human.match$new.u <- gsub("_HUMAN", "", human.match$new.u)

rat.match$new.u <- rat.match$Protein
rat.match$new.u <- gsub("..\\|.*\\|", "", rat.match$new.u)
rat.match$new.u <- gsub("_RAT", "", rat.match$new.u)
rat.match$u.match <- ifelse(rat.match$new.u %in% human.match$new.u, 
                            1, 0)
rat.u.match <- rat.match[!(rat.match$u.match == "1"),]
rat.u.match <- rat.u.match[!duplicated(rat.u.match$new.u),]


rat.h.match$new.u <- rat.h.match$Protein
rat.h.match$new.u <- gsub("..\\|.*\\|", "", rat.h.match$new.u)
rat.h.match$new.u <- gsub("_RAT", "", rat.h.match$new.u)
rat.h.match$u.match <- ifelse(rat.h.match$new.u %in% human.match$new.u, 
                              1, 0)
rat.h.u.match <- rat.h.match[!(rat.h.match$u.match == "1"),]
rat.h.u.match <- rat.h.u.match[!duplicated(rat.h.u.match$new.u),]


r <- data.frame("Peptide" = rat.u.match[1:5,2])
r$Intensity <- rat.u.match[1:5,12]
r$group <- "rat"

rh <- data.frame("Peptide" = rat.h.u.match[1:5,2])
rh$Intensity <- rat.h.u.match[1:5,12]
rh$group <- "rat.h"
r.rh <- rbind(r, rh)

ggplot(r.rh, aes(x = Peptide, y = Intensity, fill = group)) +
  geom_bar(stat="identity", color = "black", position=position_dodge()) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_classic() +
  scale_y_continuous(expand = c(0,0)) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1, size = 20),
        axis.text.y = element_text(hjust = 1, size = 20)) +
  guides(fill = F) +
  scale_fill_manual(values=c('#008000','#FF0080'))
