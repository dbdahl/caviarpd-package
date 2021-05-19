iris.dis <- dist(iris[,-5])
# In practice, use at least 100 samples and multiple cores. Less here for fast-running examples.
iris.masses <- select.masses(iris.dis, ncl.range=c(3,6), nSamplesSearch=10, nSamples=10, nCores=1)
single.mass(masses=iris.masses, distance=iris.dis, nSamples=10, nCores=1)
single.mass(masses=seq(.5, 2, by=.25), distance=iris.dis, nSamples=10, nCores=1)
