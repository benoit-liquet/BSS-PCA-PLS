#################
#### Table for result according to noise
## 3 tables: first one frequency in main text, then MSE appendix and snsibility and specificity
################

load("result_SIMU1_sample_size100.Rdata")

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.nois100 <- cbind(res.bss,res.spls)
result$SNR

load("result_SIMU1_sample_size200.Rdata")

res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.nois200 <- cbind(res.bss,res.spls)
result$SNR

load("result_SIMU1_sample_size500.Rdata")
res.spls <- cbind(result$spls,colMeans(result$MSE.spls),colMeans(result$sens.spls),colMeans(result$spe.spls))
res.bss <- cbind(result$cpls,colMeans(result$MSE.Bss),colMeans(result$sens.bss),colMeans(result$spe.bss))
res.nois500 <- cbind(res.bss,res.spls)
result$SNR


tab.Freq <- cbind(res.nois100[,c(1)],res.nois200[,c(1)],res.nois500[,c(1)],
      res.nois100[,c(5)],res.nois200[,c(5)],res.nois500[,c(5)])

tab.MSE <- cbind(res.nois100[,c(2)],res.nois200[,c(2)],res.nois500[,c(2)],
                  res.nois100[,c(6)],res.nois200[,c(6)],res.nois500[,c(6)])



tab.spe <- cbind(res.nois100[,c(3,4)],res.nois200[,c(3,4)],res.nois500[,c(3,4)],
                 res.nois100[,c(7,8)],res.nois200[,c(7,8)],res.nois500[,c(7,8)])


library(xtable)

tab.MSE <- xtable(tab.Freq,digits = 0,caption = "Number of times sparse PLS and BSS-PLS retrieve the true best subset for different model size  over 100 runs according to different noisel level")
print(tab.MSE, caption.placement = "top")

tab.MSE <- xtable(tab.MSE,digits = 2,caption = "Number of times sparse PLS and BSS-PLS retrieve the true best subset for different model size  over 100 runs according to different noisel level")
print(tab.MSE, caption.placement = "top")

tab.spe <- xtable(tab.spe,digits = 2,caption = "Sensibility (sens) and specificity (spe)  over 100 runs according to different sample size")
print(tab.spe, caption.placement = "top")
