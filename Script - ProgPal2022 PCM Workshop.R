rm(list = ls(a=T))

library(caper)
library(MCMCglmm)
library(geiger)
# library(evoldiver)


dir <- "~/Dropbox/Presentations/ProgPal 2022 - Lincoln/R workshop"

# set working directory to the folder where the data are
setwd(dir)

# read in data file as a data.frame object
df <- read.delim("data.txt", stringsAsFactors = F, header = T)
rownames(df) <- df$Taxon

# read in tree
tr <- read.nexus("tree.nex")


# Phylogenetic GLMs using caper package ----
# set up comparative data object
dt.pgls <- comparative.data(tr, df, Taxon, vcv=TRUE, vcv.dim=3)

# fit pGLS model
pgls.1 <- pgls(L_PreOrbit ~ W_Sk, dt.pgls, lambda='ML')
# print out model summary
summary(pgls.1)
# extract extract estimated lambda value from pGLS model
lambda <- pgls.1$param[2]
# rescale branches by lambda
tr.lambda <- rescale(tr, "lambda", lambda)

pdf("treeTransform.pdf")
par(mfrow = c(1,2))
# plot unscaled tree
plot(tr, show.tip.label = F)
axisPhylo()
# plot lambda-scaled tree
plot(tr.lambda, show.tip.label = F)
axisPhylo()
layout(1) # reset graphics layout
dev.off() 


# Phylogenetic GLMM using MCMCglmm package ----
# Load MCMC R package
library(MCMCglmm)

# set priors
prior1 <- list(R = list(V=1, nu=0.002), G=list(G1 = list(V = 1, nu = 1, alpha.mu = 0, alpha.V = 25^2)))
# set number of iterations
nitt <- 5E05
# set sampling interval
thin <- 1E03
# set burnin
burnin <- 1E05

# create inverse of VCV matrix
invA <- inverseA(tr, scale=F)$Ainv
# fit MCMCglmm
pglmm.1 <- MCMCglmm(L_PreOrbit ~ W_Sk, data = df, random = ~ Taxon,
         ginverse=list(Taxon=invA), family = "gaussian", prior=prior1,
         nitt=nitt, thin=thin, burnin=burnin, pl=T)
# print out model summary
summary(pglmm.1)

# save plot of the posterior distributions of the model coefficients
pdf("MCMCglmm_posteriors.pdf")
# plot MCMCglmm solutions, i.e., model coefficients
plot(pglmm.1$Sol)
dev.off()
