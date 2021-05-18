pkgname <- "caviarpd"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "caviarpd-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('caviarpd')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("caviarPD")
### * caviarPD

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: caviarPD
### Title: Cluster Analysis via Random Partition Distributions
### Aliases: caviarPD

### ** Examples

iris.dis <- dist(iris[,-5])
caviarPD(distance=iris.dis, nSamples=10)
caviarPD(distance=iris.dis, mass=0.75, loss="binder", nSamples=10, maxNClusters=3)
# In practice the user should use at least 100 samples, but for ease of testing we use less here.




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("caviarPD", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("jaccard")
### * jaccard

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: jaccard
### Title: Jaccard Distance for Categorical Attributes
### Aliases: jaccard

### ** Examples

jaccard(npk[,-5])
jaccard(warpbreaks[,-1])




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("jaccard", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("loss.indexes")
### * loss.indexes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: loss.indexes
### Title: Accuracy for a List of Clustering Estimates through Loss
###   Functions
### Aliases: loss.indexes

### ** Examples

iris.dis <- dist(iris[,-5])
iris.truth <- as.numeric(iris[,5])
est1 <- caviarPD(distance=iris.dis, mass=0.9, nSamples=1000, loss='binder')
est2 <- caviarPD(distance=iris.dis, mass=2.0, loss="VI")
loss.indexes(list(est1, est2), iris.truth)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("loss.indexes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("select.masses")
### * select.masses

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: select.masses
### Title: Mass Parameter Selection for the CaviarPD Procedure
### Aliases: select.masses

### ** Examples

tooth.dis <- dist(scale(ToothGrowth[,-2]))
# In practice, use at least 100 samples and multiple cores. Less here for fast-running examples.
select.masses(tooth.dis, ncl.range=c(2,4), nSamples=10, nSamplesFinal=10, nCores=1)
iris.dis <- dist(iris[,-5])
select.masses(iris.dis, ncl.range=c(3,6), single=TRUE, nSamples=10, nSamplesFinal=10, nCores=1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("select.masses", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("single.mass")
### * single.mass

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: single.mass
### Title: Single Mass Parameter Selection for the CaviarPD Procedure
### Aliases: single.mass

### ** Examples

iris.dis <- dist(iris[,-5])
# In practice, use at least 100 samples and multiple cores. Less here for fast-running examples.
iris.masses <- select.masses(iris.dis, ncl.range=c(3,6), nSamples=10, nSamplesFinal=10, nCores=1)
single.mass(masses=iris.masses, distance=iris.dis, nSamples=10, nCores=1)
single.mass(masses=seq(.5, 2, by=.25), distance=iris.dis, nSamples=10, nCores=1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("single.mass", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
