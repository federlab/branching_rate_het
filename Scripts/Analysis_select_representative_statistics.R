#
library(ape)
library(parallel)
library(poweRbal)

#
# list all tree balance statistics
applicable.stats = getAllTSS(n = c(50,1250),not_only_bin = TRUE,types = c("bali","imbali"))
applicable.stats = setdiff(applicable.stats,c("MeanIp","TotIp"))
balance.stats = getAllTSS(n = c(50,1250),not_only_bin = TRUE,types = "bali")
imbalance.stats = getAllTSS(n = c(50,1250),not_only_bin = TRUE,types = "imbali")

#
empirical.nulls= readRDS("Simulated_data/stats/J1_Sackin_nulls.RDS")
idx = which(empirical.nulls$meta$n==1250&
              empirical.nulls$meta$N==6250)

#
pairwise.corr = do.call(rbind,lapply(
  combn(x = c("J1",applicable.stats),m = 2,simplify = FALSE),
  function(ii){
    this.test = cor.test(do.call(c,lapply(idx,function(i){
      empirical.nulls$null[[i]][[ii[1]]]
    })),
    do.call(c,lapply(idx,function(i){
      empirical.nulls$null[[i]][[ii[2]]]
    })))
    data.frame(stat1=ii[1],stat2=ii[2],
               rho=this.test$estimate,p.value=this.test$p.value)
  }
))

#
seed.stats = c(J1="J1",Sackin="Sackin")

#
stats.cluster = lapply(seed.stats,function(x){
  df = pairwise.corr[pairwise.corr$stat1==x|pairwise.corr$stat2==x,]
  df = df[abs(df$rho)>=0.75,]
  unique(c(x,df$stat1,df$stat2))
})

remainder.stats = setdiff(applicable.stats,do.call(c,stats.cluster))
while(length(remainder.stats)>0){
  new.cluster = lapply(remainder.stats,function(x){
    df = pairwise.corr[pairwise.corr$stat1==x|pairwise.corr$stat2==x,]
    df = df[abs(df$rho)>=0.75,]
    unique(c(x,df$stat1,df$stat2))
  })
  cluster.size = sapply(new.cluster,length)
  stats.cluster[[remainder.stats[which.max(cluster.size)]]] = new.cluster[[which.max(cluster.size)]]
  remainder.stats = setdiff(applicable.stats,do.call(c,stats.cluster))
}
