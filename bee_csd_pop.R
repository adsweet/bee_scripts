library(ape)
library(reshape2)
library(matrixcalc)
library(ggplot2)
library(adegenet)
library(pegas)
library(mmod)
library(poppr)
library(apex)
library(related)

# Read in tables with population and relatedness information
bros <- read.table("brother_IDs.txt",header=T)
pops <- read.csv("bee_csd_pops.csv",header=T)

# Read in multiple sequence alignment and compute genetic distances
ex8 <- read.dna("bee_csd_ex8_1.0_k31_unique_alignment.fasta", format="fasta")
ex8.dist <- dist.dna(ex8, model = "raw", pairwise.deletion = T, as.matrix = T)
ex8.dist.tri <- lower.triangle(ex8.dist)
ex8.dist.tri[upper.tri(ex8.dist.tri, diag=T)] <- NA

ex8.dist.mat <- setNames(melt(ex8.dist.tri), c('One', 'Two', 'Distance'))
ex8.dist.mat <- na.omit(ex8.dist.mat)
write.csv(ex8.dist.mat, "bee_csd_8_alignment_dist.csv") # Write genetic distance matrix

# Population genetic statistics
## Read in alignments
ex8.dat <- read.multiFASTA("bee_csd_ex8_1.0_k31_unique_alignment.fasta")
ex7.dat <- read.multiFASTA("bee_csd_ex7_1.0_k31_unique_alignment.fasta")
csd.dat <- read.multiFASTA(c("bee_csd_ex8_1.0_k31_unique_alignment.fasta", "bee_csd_ex7_1.0_k31_unique_alignment.fasta"))
(setLocusNames(ex8.dat) <- "csd_ex8")
(setLocusNames(ex7.dat) <- "csd_ex7")
(setLocusNames(csd.dat) <- c("csd_ex8","csd_ex7"))

## Convert to Genind objects
ex8.gid <- multidna2genind(ex8.dat, mlst = T)
ex7.gid <- multidna2genind(ex7.dat, mlst = T)
csd.gid <- multidna2genind(csd.dat, mlst = T)

## Set strata
strata(ex8.gid) <- pops
setPop(ex8.gid) <- ~family
strata(csd.gid) <- pops
setPop(csd.gid) <- ~family

## Compute statistics
diff_stats(ex8.gid)
Phi_st_Meirmans(ex8.gid)
pairwise_D(ex8.gid, linearized = F, hsht_mean = "harmonic")
hap.div(ex8)

diff_stats(csd.gid)
Phi_st_Meirmans(csd.gid)
pairwise_D(csd.gid, linearized = F, hsht_mean = "harmonic")
Hs(csd.gid)
csd.dist <- dist.multidna(csd.dat, pool = T)
poppr.amova(csd.gid, ~family,nperm = 100)


## Calculate relatedness
pop.tab <- as.data.frame(csd.gid@tab)
pop.rel <- readgenotypedata(pop.tab)
hap <- haplotype(ex8)
out <- coancestry(pop.rel$gdata, lynchli = 1, lynchrd = 1, quellergt = 1, ritland = 1, wang = 1)

## Plot pairwise distances

ex8.full <- read.csv("bee_csd_8_alignment_dist_full.csv", header=T)

ggplot(ex8.full, aes(same_fam, Distance, fill=same_fam))+
  geom_boxplot()+
  geom_point(size = 1, shape = 21, position = position_jitterdodge())

## Test for differences in genetic distance within and between families (brothers)
t.test(Distance ~ same_fam, data=ex8.full)
wilcox.test(Distance ~ same_fam, data=ex8.full)
