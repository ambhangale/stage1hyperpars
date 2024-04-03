## Aditi M. Bhangale
## Last updated: 3 April 2024

# Rchaeology for runsim-sim1to3

library(doSNOW)

foo <- list(1:10, 2:11, 3:12, 4:13, 5:14)

nClus <- parallel::detectCores() - 1
cl <- makeCluster(nClus)
registerDoSNOW(cl)

bar <- foreach(list_element = 1:length(foo)) %dopar% {
                                    
        out <- mean(foo[[list_element]])
        return(out)
 }

bar

stopCluster(cl)

