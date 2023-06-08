rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)
library(jagsUI)
library(ggpubr)
library(mcmcplots)
library(reshape2)
library(grid)  


################################################################################
# Run the dilution assay model to get a value for phi

# Dilution assay data:
# expected number of molecules in dilution assay
n_mol <- c(10, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
n_mol=n_mol*1e5

# number of PCR reps run per dilution assay
n_pcr_reps <- c(8, 8, 8, 8, 8, 8, 8, 8)

# number of positive amplifications from the PCR reps
n_pos      <- c(8, 8, 8, 8, 7, 3, 2, 0)

N <- length(n_mol)

# Assay dilution JAGS model
mod.assay <- "model
  {
  for(i in 1:N) {
    n_pos[i] ~ dbin(mu[i], n_pcr_reps[i])
    mu[i] <- 1 - (1 - p[i]) ^ n[i]
    n[i] ~ dpois(n_mol[i])
    logit(p[i]) <- lp
  }
  
  lp ~ dnorm(0, 0.00001)
  logit(p_out) <- lp
}"
write(mod.assay, "JAGS_ProbAmplification_Model.txt")

out.assay <- jags(model = "JAGS_ProbAmplification_Model.txt",
                  data = list(N=N, n_mol = n_mol, n_pos = n_pos, 
                              n_pcr_reps = n_pcr_reps),
                  inits = function() list(lp = runif(1,-4.5, -3.8)),   
                  param = c("lp", "p_out"),
                  n.chains = 3,
                  n.iter =20000,
                  n.burnin = 10000)

# Print and plot results

all.sum <- out.assay$summary
all.sum[, c(1,2,3,7,8)]

# extract value for p and its standard deviation
p.mean <- all.sum[2, 1]
p.se <- all.sum[2, 2]

# extract value for logit p and its standard deviation
logit.p.mean <- all.sum[1, 1]
logit.p.se <- all.sum[1, 2]



#Plot results

xx <- seq(min(n_mol) , max(n_mol), 1)  
crve <- data.frame(y = (1 - (1 - p.mean) ^ xx), x = xx)
dd <- data.frame(n_mol = n_mol, n_pcr_reps = n_pcr_reps, n_pos = n_pos)
dd=rbind(dd, data.frame(n_mol=100, n_pcr_reps=8, n_pos=7))
dd$av_prob <- dd$n_pos / dd$n_pcr_reps 

p1=ggplot(dd, aes(y = av_prob*100, x = n_mol)) +
  geom_point(size = 5) +
  scale_x_log10() +
  xlab("\nExpected number of molecules in PCR replicate") +
  ylab("Proportion of successful amplifications\n") +
  geom_line(data = crve, aes(y = y*100, x = x), lwd = 1) +
  theme_classic(18) 

p.out <- data.frame(p.post = out.assay$sims.list[2]) 
vp <- viewport(width = 0.4, height = 0.4, x = 0.75, y = 0.4)
p2=ggplot(p.out, aes(x = p_out)) + geom_density(fill = "grey") +
  xlab(expression(italic(phi))) + ylab("Density") +
  theme_classic(14) 

print(p1)
print(p2, vp = vp)

################################################################################




################################################################################
# Concentration and Sensitivity Model

mdat2 <- fread("Input_data_SampleIDs.csv")

Latrines=mdat2[SiteType!="Non-latrine habitat" & SiteType!="STQ positive camera trap"]

sites_names <- Latrines$LocationCode
sites_fact <- as.numeric(factor(sites_names))
samples_fact <- as.numeric(factor(Latrines[, Sample]))
N_sites <- length(unique(sites_names))
N_samples <- length(Latrines[, unique(Sample)])

pcr <- array(NA, dim=c(N_samples, 2, 3))
sites <- vector("integer", length = N_samples)
samples <- vector("integer", length = N_samples)

distance_m <- samples <- vector("integer", length = N_samples)
siteType <- samples <- vector("integer", length = N_samples)
for(s in seq_len(N_samples)) { # sample
  for(j in 1:2) { # qPCR
    for(z in 1:3) { # extracts
      pcr[s,j,z] <- Latrines[(s-1)*3+z, ][[5+j]]
      if(j==1 & z==1) {
        samples[s] <- samples_fact[(s-1)*3+z]
        sites[s] <- sites_fact[(s-1)*3+z]
        distance_m[s] <- Latrines$DistanceFromLatrine[(s-1)*3+z]
        siteType[s] <- Latrines$SiteType[(s-1)*3+z]
      }
    }
  }
}

N_pcr_reps <- 2 #number of PCR replicates per extract
N_extracts <- 3 #number of extracts per sample


sample_vol=matrix(rep(250, (N_samples * N_extracts)), nrow = N_samples) #actually sample weight in this case (mg)
pcr_dilution <- 2/100


# Prior for prob of amplification of one DNA mol estimated with dilution series
logit.p.mean <- -3.422643 
logit.p.se <- 0.3730191

setkey(Latrines, LocationCode)
site.dt <- Latrines[Latrines[, unique(LocationCode)], mult="first", .(LocationCode, SiteType)]
act_sites <- which(site.dt[, SiteType] == "Active Latrine")
inact_sites <- which(site.dt[, SiteType] == "Inactive Latrine")

#Can load model results to avoid having to rerun it
load("mod.Rda")


#Calculate sensitivity for greater number of samples, 
# Run the model in JAGS with six chains
mod <- jags(model="./JAGS_Sensitivity_Model.txt",
                    data=list(y=pcr, d=distance_m, sites=sites,
                              N_bags=N_samples, N_pcr_reps=N_pcr_reps,
                              N_extracts=N_extracts,
                              N_sites=N_sites,
                              logit.p.mean=logit.p.mean, logit.p.se=logit.p.se,
                              pcr_dilution=pcr_dilution, sample_vol=sample_vol,
                              N_samples=10, 
                              act_sites=act_sites, inact_sites=inact_sites,
                              N_site_types=2
                    ),
                    inits=function() list(n_mol_pcr=pcr, 
                                          beta_d=-abs(rnorm(1,0,0.1)),
                                          eps=rnorm(N_samples,0,1)
                                          ),
                    param=c("concentration", "p_out", "r", "beta_d", "sd_eps",
                            "eps", 
                            "mu", "n_mol_sample", "n_mol_pcr", 
                            "sens", "theta", "eta"),
                    n.chains=6,
                    n.iter = 350000,
                    n.thin = 300,
                    n.burnin=20000,
            parallel = TRUE)
save(mod, file="mod.Rda")


mod.sum <- mod$summary


#Extract and plot concentration results
site.dt[, Concentration:= mod.sum[grepl("^concen", row.names(mod.sum)), "mean"]] #concentration at latrine
site.dt[, Concentration.llim:= mod.sum[grepl("^concen", row.names(mod.sum)), "2.5%"]] #concentration at latrine
site.dt[, Concentration.ulim:= mod.sum[grepl("^concen", row.names(mod.sum)), "97.5%"]] #concentration at latrine

site.dt[, Concentration.1:=Concentration + mod.sum[grepl("beta_d", row.names(mod.sum)), "mean"]] #concentration at 1m
site.dt[, Concentration.1.llim:=Concentration.llim + mod.sum[grepl("beta_d", row.names(mod.sum)), "2.5%"]] #concentration at 1m
site.dt[, Concentration.1.ulim:=Concentration.ulim + mod.sum[grepl("beta_d", row.names(mod.sum)), "97.5%"]] #concentration at 1m

site.dt[, .(mean0=mean(exp(Concentration)), llim0=mean(exp(Concentration.llim)), ulim0=mean(exp(Concentration.ulim))), by=SiteType] #average across site types
site.dt[, .(mean1=mean(exp(Concentration.1)), llim1=mean(exp(Concentration.1.llim)), ulim1=mean(exp(Concentration.1.ulim))), by=SiteType] #average across site types
site.dt[, .(mean0=mean(exp(Concentration)), llim0=mean(exp(Concentration.llim)), ulim0=mean(exp(Concentration.ulim))), by=LocationCode] #average by site
site.dt[, .(mean1=mean(exp(Concentration.1)), llim1=mean(exp(Concentration.1.llim)), ulim1=mean(exp(Concentration.1.ulim))), by=LocationCode] #average by site

site.dt[, Concentration.e:=exp(Concentration)]
site.dt[, Concentration.1.e:=exp(Concentration.1)]

p1=ggplot(site.dt, aes(SiteType, Concentration.e)) + 
  geom_boxplot(outlier.shape = NA, fill="forestgreen") + geom_jitter(width=0.05, size=2) +
  theme_classic() + annotate("text", label="At latrine", x=0.65, y=80, size=6) +
  scale_x_discrete(labels=c("Active latrines", "Inactive latrines")) +
  stat_summary(fun=mean, geom='point', shape=18, size=8) +
  ylab("Concentration (molecules/mg)") + xlab("") +
  theme(text=element_text(size=15))

p2=ggplot(site.dt, aes(SiteType, Concentration.1.e)) + 
  geom_boxplot(outlier.shape = NA, fill="forestgreen") + geom_jitter(width=0.05, size=2) +
  theme_classic() + annotate("text", label="1m from\nlatrine", x=0.65, y=1.3, size=6) +
  scale_x_discrete(labels=c("Active latrines", "Inactive latrines")) +
  stat_summary(fun=mean, geom='point', shape=18, size=8) +
  ylab("Concentration (molecules/mg)") + xlab("") +
  theme(text=element_text(size=15))

ggarrange(p1, p2, nrow=2)



#Beta
beta=(mod.sum["beta_d",])
beta
100*(exp(beta["mean"])-1) #percent change in DNA concentration at 1m from latrine site
100*(exp(beta["2.5%"])-1)
100*(exp(beta["97.5%"])-1)


#Extract Dispersion Results
r=mod.sum["r",]
r


#Extract and plot Sensitivity Results
sens <- as.data.frame(mod.sum[grep(pattern = "sens", x = row.names(mod.sum)),])
sens$Distance <- rep(0:1, each=20)
sens$nSamples <- rep(1:10, each=2)
sens$Type <- rep(c("Active", "Inactive"), 20)
names(sens)[c(3, 5, 7)] <- c("llim", "Median", "ulim")

lat.labs <- c("Active latrines", "Inactive latrines")
names(lat.labs) <- c("Active", "Inactive")
dist.labs <- c("At latrine", "1m from latrine")
names(dist.labs) <- c("0", "1")

ggplot(sens, aes(nSamples, Median)) + 
  geom_line(lwd=1) +
  geom_hline(yintercept=0.95, lty=2) +
  geom_ribbon(aes(ymin=llim, ymax=ulim), alpha=0.3, fill="red") +
  geom_ribbon(aes(ymin=`25%`, ymax=`75%`), alpha=0.3, fill="red") +
  facet_grid(Type~Distance, labeller = labeller(Type = lat.labs, Distance = dist.labs)) + 
  ylab("Probability of Detection") +
  xlab("Number of samples analysed") +
  scale_x_continuous(breaks=c(2,4,6,8,10)) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
        panel.background = element_rect(colour = "black", linewidth=1),
        axis.line = element_line(colour="black"),
        panel.grid.major = element_line(colour="#f0f0f0"),
        text = element_text(size=15))

