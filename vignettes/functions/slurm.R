
envir <- commandArgs(trailingOnly = T)[1]
i <- as.integer(commandArgs(trailingOnly = T)[2])

load(envir)

.libPaths(paths)

toattach <- s.pckgs

for (p in toattach ) library(p, character.only = T)

out <- do.call(FUN, c(unname(X[i]),args ))

save(out, file = sprintf("%s-%d",envir,i))
