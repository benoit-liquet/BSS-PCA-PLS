load("Result_PCA_n50.RData")
res <- result$selection
freq <- result$freq
load("Result_PCA_n100.RData")
res <- rbind(res,result$selection)
freq <- cbind(freq,result$freq)
load("Result_PCA_n300.RData")
res <- rbind(res,result$selection)
freq <- cbind(freq,result$freq)

library(xtable)
rownames(res) <- rep(c("PCA","BSS-PCA","sPCA","SPCA"),3)
tab.res <- xtable(res,digits = c(0,2,2,2,2,2,2),caption = "Number of times sparse PLS and BSS-PLS retrieve the true best subset for different model size  over 100 runs according to different noisel level")
print(tab.res, caption.placement = "top")


tab.freq <- xtable(freq,digits = c(0),caption = "Number of times sparse PLS and BSS-PLS retrieve the true best subset for different model size  over 100 runs according to different noisel level")
print(tab.freq, caption.placement = "top")




