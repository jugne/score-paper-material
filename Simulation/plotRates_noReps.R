######################################################
######################################################
# Here the inferred mean coalescent and migration
# rate ratios are plotted
######################################################
######################################################
library(ggplot2)
library(coda)


# clear workspace
rm(list = ls())



# Set the directory to the directory of the file
this.dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(this.dir)

# read in the true rates
true_rates_file <- "rates.csv"
true.rates <- read.table(true_rates_file, header=TRUE, sep=",")

# use the matlab standard colors to plot
col0 <- rgb(red=0.0, green=0.4470,blue=0.7410)
col1 <- rgb(red=0.8500, green=0.3250,blue=0.0980)
col2 <- rgb(red=0.9290, green=0.6940,blue=0.1250)
col4 <- rgb(red=0.4660, green=0.6740,blue=0.1880)
col3 <- rgb(red=0.3010, green=0.7450,blue=0.9330)


gen_information = c("low", "mixed", "high")


for (g in seq(1,length(gen_information))){
  print(gen_information[[g]])
# get the names of all output files of the first replicate
  log <- list.files(path="./inference/", pattern=paste("inf_", gen_information[[g]], ".*.log", sep=""), full.names = TRUE)
first = T
c = 0;
for (i in seq(1,length(log))){
  # read in the log file
  t <- read.table(log[[i]], header=TRUE, sep="\t")
  
  # take a 10% burnin
  t <- t[-seq(1,ceiling(length(t$posterior)/5)), ]
  
  # calculate ess values
  ess <- effectiveSize(t)
  #
  post_ess = as.numeric(ess["posterior"])
  # print(post_ess)

  if (post_ess>100){

    # find the correct run number
    tmp = strsplit(strsplit(log[[i]], '_')[[1]][[3]], '\\.')[[1]][[1]]
    used.rates = true.rates[which(true.rates$run==as.numeric(tmp)),];
    
    # combine with the other replicates
    new.rea_1 = data.frame(true=used.rates$reassortment_1, estimated=mean(t$reassortmentRate.0), 
                           upper=quantile(t$reassortmentRate.0,0.975), lower=quantile(t$reassortmentRate.0,0.025) )
    new.rea_2 = data.frame(true=used.rates$reassortment_2, estimated=mean(t$reassortmentRate.1), 
                           upper=quantile(t$reassortmentRate.1,0.975), lower=quantile(t$reassortmentRate.1,0.025) )
    new.mig_1 = data.frame(true=used.rates$migration_1, estimated=mean(t$b_migrationRate.0_to_1), 
                           upper=quantile(t$b_migrationRate.0_to_1,0.975), lower=quantile(t$b_migrationRate.0_to_1,0.025) )
    new.mig_2 = data.frame(true=used.rates$migration_2, estimated=mean(t$b_migrationRate.1_to_0), 
                           upper=quantile(t$b_migrationRate.1_to_0,0.975), lower=quantile(t$b_migrationRate.1_to_0,0.025) )
    new.Ne_1 = data.frame(true=used.rates$Ne_1, estimated=mean(t$popSize.0), 
                          upper=quantile(t$popSize.0,0.975), lower=quantile(t$popSize.0,0.025))
    new.Ne_2 = data.frame(true=used.rates$Ne_2, estimated=mean(t$popSize.1), 
                          upper=quantile(t$popSize.1,0.975), lower=quantile(t$popSize.1,0.025))
    
    
    hpd.rea_1 = HPDinterval(as.mcmc(t$reassortmentRate.0))
    hpd.rea_2 = HPDinterval(as.mcmc(t$reassortmentRate.1))
    hpd.mig_1 = HPDinterval(as.mcmc(t$b_migrationRate.0_to_1))
    hpd.mig_2 = HPDinterval(as.mcmc(t$b_migrationRate.1_to_0))
    hpd.Ne_1 = HPDinterval(as.mcmc(t$popSize.0))
    hpd.Ne_2 = HPDinterval(as.mcmc(t$popSize.1))
    
    new.test.rea_1 = c(as.numeric(used.rates$reassortment_1 >= hpd.rea_1[1] & used.rates$reassortment_1 <= hpd.rea_1[2]))
    new.test.rea_2 = c(as.numeric(used.rates$reassortment_2 >= hpd.rea_2[1] & used.rates$reassortment_2 <= hpd.rea_2[2]))
    
    new.test.mig_1 = c(as.numeric(used.rates$migration_1 >= hpd.mig_1[1] & used.rates$migration_1 <= hpd.mig_1[2]))
    new.test.mig_2 = c(as.numeric(used.rates$migration_2 >= hpd.mig_2[1] & used.rates$migration_2 <= hpd.mig_2[2]))
    
    new.test.Ne_1 = c(as.numeric(used.rates$Ne_1 >= hpd.Ne_1[1] & used.rates$Ne_1 <= hpd.Ne_1[2]))
    new.test.Ne_2 = c(as.numeric(used.rates$Ne_2 >= hpd.Ne_2[1] & used.rates$Ne_2 <= hpd.Ne_2[2]))
    
    
    if (first){
      rea_1 = new.rea_1
      mig_1 = new.mig_1
      mig_2 = new.mig_2
      rea_2 = new.rea_2
      Ne_1 = new.Ne_1
      Ne_2 = new.Ne_2
      
      test.rea_1 = new.test.rea_1
      test.rea_2 = new.test.rea_2
      test.mig_1 = new.test.mig_1
      test.mig_2 = new.test.mig_2
      test.Ne_1 = new.test.Ne_1
      test.Ne_2 = new.test.Ne_2
      
      
      first = F
    }else{
      rea_1 = rbind(rea_1, new.rea_1)
      rea_2 = rbind(rea_2, new.rea_2)
      mig_1 = rbind(mig_1, new.mig_1)
      mig_2 = rbind(mig_2, new.mig_2)
      Ne_1 = rbind(Ne_1, new.Ne_1)
      Ne_2 = rbind(Ne_2, new.Ne_2)
      
      test.rea_1 = rbind(test.rea_1, new.test.rea_1)
      test.rea_2 = rbind(test.rea_2, new.test.rea_2)
      
      test.mig_1 = rbind(test.mig_1, new.test.mig_1)
      test.mig_2 = rbind(test.mig_2, new.test.mig_2)
      
      test.Ne_1 = rbind(test.Ne_1, new.test.Ne_1)
      test.Ne_2 = rbind(test.Ne_2, new.test.Ne_2)
      
    }
  }else{
    c = c+1
    print(i)
  }
}

y.min_rea = min(rbind(rea_1$lower, rea_2$lower))
y.max_rea = max(rbind(rea_1$upper, rea_2$upper))
y.min_Ne = min(rbind(Ne_1$lower, Ne_2$lower))
y.max_Ne = max(rbind(Ne_1$upper, Ne_2$upper))
y.min_mig = min(rbind(mig_1$lower, mig_2$lower))
y.max_mig = max(rbind(mig_1$upper, mig_2$upper))

x.min_rea = min(rbind(rea_1$true, rea_2$true))
x.max_rea = max(rbind(rea_1$true, rea_2$true))
x.min_Ne = min(rbind(Ne_1$true, Ne_2$true))
x.max_Ne = max(rbind(Ne_1$true, Ne_2$true))
x.min_mig = min(rbind(mig_1$true, mig_2$true))
x.max_mig = max(rbind(mig_1$true, mig_2$true))

p.Ne_1 <- ggplot(Ne_1)+
  geom_abline(intercept = 0, color="red")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10(limits=c(y.min_Ne, y.max_Ne)) +
  scale_x_log10(limits=c(x.min_Ne, x.max_Ne)) 
p.Ne_2 <- ggplot(Ne_2)+
  geom_abline(intercept = 0, color="red")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10(limits=c(y.min_Ne, y.max_Ne)) +
  scale_x_log10(limits=c(x.min_Ne, x.max_Ne)) 
p.rea_1 <- ggplot(rea_1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10(limits=c(y.min_rea, y.max_rea)) +
  scale_x_log10(limits=c(x.min_rea, x.max_rea)) 
p.rea_2 <- ggplot(rea_2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10(limits=c(y.min_rea, y.max_rea)) +
  scale_x_log10(limits=c(x.min_rea, x.max_rea)) 
p.mig_1 <- ggplot(mig_1)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10(limits=c(y.min_mig, y.max_mig)) +
  scale_x_log10(limits=c(x.min_mig, x.max_mig)) 
p.mig_2 <- ggplot(mig_2)+
  geom_abline(intercept = 0, color="red", linetype="dashed")+
  geom_errorbar(aes(x=true, ymin=lower, ymax=upper), colour="grey", width=0.01) +
  geom_point(aes(x=true, y=estimated), size=2) + 
  theme_minimal() +
  scale_y_log10(limits=c(y.min_mig, y.max_mig)) +
  scale_x_log10(limits=c(x.min_mig, x.max_mig)) 

# scale_y_log10(limits=c(0.01,1.0)) +
# scale_x_log10(limits=c(0.01,1.0)) 

plot(p.Ne_1)
plot(p.rea_1)
plot(p.Ne_2)
plot(p.rea_2)
plot(p.mig_1)
plot(p.mig_2)

ggsave(plot=p.Ne_1,paste("figures_logN_mig_prior/type_1_Ne_sim_",gen_information[[g]], ".pdf", sep=""),width=6, height=5)
ggsave(plot=p.rea_1,paste("figures_logN_mig_prior/type_1_rho_sim_",gen_information[[g]], ".pdf", sep=""),width=6, height=5)
ggsave(plot=p.mig_1,paste("figures_logN_mig_prior/type_1_mig_sim_",gen_information[[g]], ".pdf", sep=""),width=6, height=5)
ggsave(plot=p.Ne_2,paste("figures_logN_mig_prior/type_2_Ne_sim_",gen_information[[g]], ".pdf", sep=""),width=6, height=5)
ggsave(plot=p.rea_2,paste("figures_logN_mig_prior/type_2_rho_sim_",gen_information[[g]], ".pdf", sep=""),width=6, height=5)
ggsave(plot=p.mig_2,paste("figures_logN_mig_prior/type_2_mig_sim_",gen_information[[g]], ".pdf", sep=""),width=6, height=5)

hpd.test = data.frame("Parameter"=c('rea_1', 'rea_2', 'mig_1', 'mig_2', 'Ne_1', 'Ne_2'), 
                      "HPD coverage"=c(mean(test.rea_1), mean(test.rea_2), mean(test.mig_1), 
                                       mean(test.mig_2), mean(test.Ne_1), mean(test.Ne_2)))


write.csv(hpd.test, file = paste("HPD_test_",gen_information[[g]],".csv", sep="" ))

}


