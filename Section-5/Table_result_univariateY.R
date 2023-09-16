#################
## 1 table appendix:  MSE  and snsibility and specificity
################

load("result_SIMU_Lasso_sparsity_n_p_sparisty_SNR400401010.Rdata")
res <- cbind(result$MSE.bss,result$sens.bss,result$spe.bss,result$MSE.spls,result$sens.spls,result$spe.spls,result$MSE.lasso,result$sens.lasso,result$spe.lasso)
result1 <- c(400,40,10,10,colMeans(res))
result$SNR

load("result_SIMU_Lasso_sparsity_n_p_sparisty_SNR400401015.Rdata")
res <- cbind(result$MSE.bss,result$sens.bss,result$spe.bss,result$MSE.spls,result$sens.spls,result$spe.spls,result$MSE.lasso,result$sens.lasso,result$spe.lasso)
result2 <- c(400,40,10,15,colMeans(res))
result$SNR

load("result_SIMU_Lasso_sparsity_n_p_sparisty_SNR400403010.Rdata")
res <- cbind(result$MSE.bss,result$sens.bss,result$spe.bss,result$MSE.spls,result$sens.spls,result$spe.spls,result$MSE.lasso,result$sens.lasso,result$spe.lasso)
result3 <- c(400,40,30,10,colMeans(res))
result$SNR


load("result_SIMU_Lasso_sparsity_n_p_sparisty_SNR400403015.Rdata")
res <- cbind(result$MSE.bss,result$sens.bss,result$spe.bss,result$MSE.spls,result$sens.spls,result$spe.spls,result$MSE.lasso,result$sens.lasso,result$spe.lasso)
result4 <- c(400,40,30,15,colMeans(res))
result$SNR

load("result_SIMU_Lasso_sparsity_n_p_sparisty_SNR_new40802022.Rdata")
res <- cbind(result$MSE.bss,result$sens.bss,result$spe.bss,result$MSE.spls,result$sens.spls,result$spe.spls,result$MSE.lasso,result$sens.lasso,result$spe.lasso)
result5 <- c(40,80,20,10,colMeans(res))
result$SNR

load("result_SIMU_Lasso_sparsity_n_p_sparisty_SNR_new40802038.Rdata")
res <- cbind(result$MSE.bss,result$sens.bss,result$spe.bss,result$MSE.spls,result$sens.spls,result$spe.spls,result$MSE.lasso,result$sens.lasso,result$spe.lasso)
result6 <- c(40,80,20,15,colMeans(res))
result$SNR

load("result_SIMU_Lasso_sparsity_n_p_sparisty_SNR_new40804022.Rdata")
res <- cbind(result$MSE.bss,result$sens.bss,result$spe.bss,result$MSE.spls,result$sens.spls,result$spe.spls,result$MSE.lasso,result$sens.lasso,result$spe.lasso)
result7 <- c(40,80,40,10,colMeans(res))
result$SNR

load("result_SIMU_Lasso_sparsity_n_p_sparisty_SNR_new40804038.Rdata")
res <- cbind(result$MSE.bss,result$sens.bss,result$spe.bss,result$MSE.spls,result$sens.spls,result$spe.spls,result$MSE.lasso,result$sens.lasso,result$spe.lasso)
result8 <- c(40,80,40,20,colMeans(res))
result$SNR

tab <- rbind(result1,result2,result3,result4,result5,result6,result7,result8)
tab


tab.MSE <- xtable(tab,digits = c(0,0,0,0,0,2,2,2,2,2,2,2,2,2))
print(tab.MSE, caption.placement = "top",include.rownames=F)





