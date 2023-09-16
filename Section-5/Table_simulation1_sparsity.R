#################
#### Table for result according to noise
## 3 tables: first one frequency in main text, then MSE appendix and snsibility and specificity
################

load("result_SIMU1_sparsity3.Rdata")

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.3 <- cbind(res.bss,res.spls)

load("result_SIMU1_sparsity7.Rdata")

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.7 <- cbind(res.bss,res.spls)

load("result_SIMU1_sparsity9.Rdata")
res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.9 <- cbind(res.bss,res.spls)

load("result_SIMU1_sparsity11.Rdata")
res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.11 <- cbind(res.bss,res.spls)


tab.Freq <- cbind(res.3[,c(1)],res.7[,c(1)],res.9[,c(1)],res.11[,c(1)],
                  res.3[,c(5)],res.7[,c(5)],res.9[,c(5)],res.11[,c(5)])

tab.MSE <- cbind(res.3[,c(2)],res.7[,c(2)],res.9[,c(2)],res.11[,c(2)],
                   res.3[,c(6)],res.7[,c(6)],res.9[,c(6)],res.11[,c(6)])

tab.spe <- cbind(res.3[,c(3,4)],res.7[,c(3,4)],res.9[,c(3,4)],res.9[,c(3,4)],
                 res.3[,c(7,8)],res.7[,c(7,8)],res.9[,c(7,8)],res.11[,c(7,8)])


library(xtable)

tab.Freq <- xtable(tab.Freq,digits = 0,caption = "Number of times sparse PLS and BSS-PLS retrieve the true best subset for different model size  over 100 runs according to different noisel level")
print(tab.Freq, caption.placement = "top")

tab.MSE <- xtable(tab.MSE,digits = 2,caption = "Number of times sparse PLS and BSS-PLS retrieve the true best subset for different model size  over 100 runs according to different noisel level")
print(tab.MSE, caption.placement = "top")

tab.spe <- xtable(tab.spe,digits = 2,caption = "Sensibility (sens) and specificity (spe)  over 100 runs according to different sample size")
print(tab.spe, caption.placement = "top")
