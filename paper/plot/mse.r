#! /usr/bin/env Rscript
adaptive      = read.table("mse_adaptive.dat")
non_adaptive  = read.table("mse_noadaptive.dat")
no_importance = read.table("mse_noimportance.dat")

a = 0.0001

bin <- function(x, data){
    data_hist <- margin.table(table((round(data[[1]]/ x) * x)),1)
    return(data_hist)
}

write.table(bin(a,adaptive), file="mse_adaptive_hist.dat", quote=FALSE)
write.table(bin(a,non_adaptive), file="mse_noadaptive_hist.dat", quote=FALSE)
write.table(bin(a,no_importance), file="mse_noimportance_hist.dat", quote=FALSE)
