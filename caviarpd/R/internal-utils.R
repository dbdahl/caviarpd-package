# This is where we put small, helpful functions that are not exported.

mkSeed <- function() as.raw(sample.int(256L,16L)-1L)
