##### difference in mean m_e between chromosomes

rm(list=ls())

library(ggplot2)
library(scales)
library(naniar)
library(caTools)

# read in gIMble results
dat <- read.table(file="../Downloads/brenthis_ino_daphne.autosomes_plus_Z.gridsearch_windows_kmax_2_grid_v1_0.best.bed", 
                  header=T, sep = "\t")

# chromosomes from which to sample from
ino_chromosomes <- c("brenthis_ino.SP_BI_364.chromosome_1", "brenthis_ino.SP_BI_364.chromosome_2", 
                   "brenthis_ino.SP_BI_364.chromosome_3", "brenthis_ino.SP_BI_364.chromosome_4",
                   "brenthis_ino.SP_BI_364.chromosome_5", "brenthis_ino.SP_BI_364.chromosome_6", 
                   "brenthis_ino.SP_BI_364.chromosome_7", "brenthis_ino.SP_BI_364.chromosome_8", 
                   "brenthis_ino.SP_BI_364.chromosome_9", "brenthis_ino.SP_BI_364.chromosome_10",
                   "brenthis_ino.SP_BI_364.chromosome_11", "brenthis_ino.SP_BI_364.chromosome_12", 
                   "brenthis_ino.SP_BI_364.chromosome_13", "brenthis_ino.SP_BI_364.chromosome_14")

null_values <- c()
for (i in seq(1, 100000)){ # do enough permutations to get all 3003 possible labelings
r_chromosomes <- sample(ino_chromosomes, 6)
dat$rearranged <- "0"
dat$rearranged[dat$sequence %in% r_chromosomes] <- "1"
dat$rearranged <- as.factor(dat$rearranged)
replicate <- mean(subset(dat, rearranged =="0")$me) - mean(subset(dat, rearranged =="1")$me)
null_values <- append(null_values, replicate)
}

null_values <- unique(null_values)
length(null_values) # should be 3003

hist(null_values)
mean(null_values) # should be very close to zero

null_values <- sort(null_values, decreasing = TRUE)

# observed value is 1.082e-7

null_values[floor(0.0005 * 3003)]
null_values[floor(0.001 * 3003)]
null_values[floor(0.005 * 3003)] ### this generates a smaller test statistic
null_values[floor(0.01 * 3003)]
null_values[floor(0.05 * 3003)]

#################################################################################
### difference in mean m_e between windows that are close or far from rearrangement points


rm(list=ls())

library(ggplot2)
library(scales)
library(naniar)
library(caTools)

dat <- read.table(file="../Downloads/brenthis_ino_daphne.autosomes_plus_Z.gridsearch_windows_kmax_2_grid_v1_0.best.bed", 
                  header=T, sep = "\t")

# set these, which remain constant
r_chromosomes <- c("brenthis_ino.SP_BI_364.chromosome_1", "brenthis_ino.SP_BI_364.chromosome_3", 
                   "brenthis_ino.SP_BI_364.chromosome_5", "brenthis_ino.SP_BI_364.chromosome_8",
                   "brenthis_ino.SP_BI_364.chromosome_9", "brenthis_ino.SP_BI_364.chromosome_14")

# mark chromosomes as rearranged
dat$rearranged <- "0"
dat$rearranged[dat$sequence %in% r_chromosomes] <- "1"
dat$rearranged <- as.factor(dat$rearranged)

# subset data and add an ID column
temp_dat <- subset(dat, rearranged == "1")
temp_dat$ID <- seq.int(nrow(temp_dat))

# a reminder of what the real data looks like
#tightly_linked <- c(25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 646, 
#                    647, 648, 649, 650, 651, 652, 653, 654, 655, 656, 657, 658, 
#                    659, 660, 710, 711, 712, 2667, 2668, 2669, 2670, 2671, 
#                    2672, 1773, 1774, 1775, 1776, 1777, 1778, 1779, 1780, 1781, 
#                    1782, 1783, 1784, 1803, 1804, 1832, 1833, 1834, 1835, 1836, 
#                    1837, 1838, 1839, 1840, 1841, 1842, 1843, 1974, 1975, 1976, 
#                    1977, 1978, 1979, 1980, 1981, 1168, 1169, 1170, 1171, 1172)

# a function that allow adjacent windows to roll over
roll_over <- function(x) {
  if (x > 1249){
    return (x - 1249)
  } else {
    return (x)
  }
}

more_null_values <- c()
sampled <- c()
for (i in seq(1, 100000)){ # do as many samples to get a good approximation
  simmed_tl <- c()
  for (j in c(12, 14, 2, 5, 11, 1, 11, 7, 4)){ # corresponds to sizes in real data
    sampled_i <- sample(temp_dat$ID, 1)
    sampled_range <- seq(sampled_i, sampled_i + j)
    sampled_range <- sapply(sampled_range, roll_over)
    simmed_tl <- append(simmed_tl, sampled_range)
    }
  temp_dat$tl <- "0"
  temp_dat$tl[temp_dat$ID %in% simmed_tl] <- "1"
  temp_dat$tl <- as.factor(temp_dat$tl)
  replicate <- mean(subset(temp_dat, rearranged =="1" & tl=="0")$me) - mean(subset(temp_dat, rearranged =="1" & tl=="1")$me)
  hits <- nrow(subset(temp_dat, rearranged =="1" & tl=="1"))
  if (hits == 76){ # do this to condition on non-overlapping windows
    more_null_values <- append(more_null_values, replicate)
    sampled <- append(sampled, simmed_tl)
  }
}

length(more_null_values) # will be < 100000
mean(more_null_values) # should be close to zero
hist(more_null_values, breaks=50)
hist(sampled, breaks=125) # chromosome ends should have been well sampled

more_null_values <- sort(more_null_values, decreasing = TRUE)

# observed value is 7.662e-8 

lmnv <- length(more_null_values)
more_null_values[floor(lmnv * 0.0001)]
more_null_values[floor(lmnv * 0.0005)] ### this generates a smaller test statistic
more_null_values[floor(lmnv * 0.001)]
more_null_values[floor(lmnv * 0.005)]
more_null_values[floor(lmnv * 0.01)]
more_null_values[floor(lmnv * 0.05)]

#################################################################################
### difference in barrier frequency between chromosomes

rm(list=ls())

library(ggplot2)
library(scales)
library(naniar)
library(caTools)

dat_me <- read.table(file="../Downloads/brenthis_ino_daphne.autosomes_plus_Z.gridsearch_windows_kmax_2_grid_v1_0.fixed_param.me.bed", 
                     header=T, sep = "\t")

ino_chromosomes <- c("brenthis_ino.SP_BI_364.chromosome_1", "brenthis_ino.SP_BI_364.chromosome_2", 
                     "brenthis_ino.SP_BI_364.chromosome_3", "brenthis_ino.SP_BI_364.chromosome_4",
                     "brenthis_ino.SP_BI_364.chromosome_5", "brenthis_ino.SP_BI_364.chromosome_6", 
                     "brenthis_ino.SP_BI_364.chromosome_7", "brenthis_ino.SP_BI_364.chromosome_8", 
                     "brenthis_ino.SP_BI_364.chromosome_9", "brenthis_ino.SP_BI_364.chromosome_10",
                     "brenthis_ino.SP_BI_364.chromosome_11", "brenthis_ino.SP_BI_364.chromosome_12", 
                     "brenthis_ino.SP_BI_364.chromosome_13", "brenthis_ino.SP_BI_364.chromosome_14")

# calculate deltaB and a binary measure of whether a barrier or not
dat_me$deltaB <- dat_me$me_0.0 - dat_me$me_1.75e.07
dat_me$binaryB <- 0
dat_me$binaryB[dat_me$deltaB > 0] <- 1

null_values <- c()
for (i in seq(1, 50000)){ # do enough permutations to sample each labeling
  r_chromosomes <- sample(ino_chromosomes, 6)
  dat_me$rearranged <- "0"
  dat_me$rearranged[dat_me$sequence %in% r_chromosomes] <- "1"
  dat_me$rearranged <- as.factor(dat_me$rearranged)
  replicate <- mean(subset(dat_me, rearranged =="1")$binaryB) - mean(subset(dat_me, rearranged =="0")$binaryB)
  r_chromosomes <- append(r_chromosomes, replicate)
  r_chromosomes <- list(sort(r_chromosomes))
  null_values <- append(null_values, r_chromosomes)
}

# I have recorded both the labeling and the result
# this is because some labelings give the same result
# so I cannot just use unique on the results

null_values <- unique(null_values)
length(null_values) # should be 3003
null_values <- sapply(null_values, "[[",1) # get rid of labelings
null_values <- as.double(null_values)

hist(null_values)
mean(null_values)

null_values <- sort(null_values, decreasing=TRUE)

# observed value is 0.1764

lnv <- length(null_values)
null_values[floor(0.0005 * lnv)]
null_values[floor(0.001 * lnv)]
null_values[floor(0.005 * lnv)]
null_values[floor(0.01 * lnv)] ### this generates a smaller test statistic
null_values[floor(0.05 * lnv)]
