tooth.dis <- dist(scale(ToothGrowth[,-2]))
# In practice, use at least 100 samples and multiple cores. Less here for fast-running examples.
select.masses(tooth.dis, ncl.range=c(2,4), nSamplesSearch=10, nSamples=10, nCores=1)
iris.dis <- dist(iris[,-5])
select.masses(iris.dis, ncl.range=c(3,6), single=TRUE, nSamplesSearch=10, nSamples=10, nCores=1)
