################################################################################
#
# Weight Trial. 
# 
# Author: Alan Vazquez
# Affiliation: Tecnologico de Monterrey
# Email: alanrvazquez@tec.mx
#
################################################################################

library(dplyr)
library(reshape2)
library(tidyr)
library(gurobi) # Load Gurobi solver.
source("functions/plot_theme.R") # To enhance plots.

# Read data.
trial.data = read.csv("data/DATA.csv")
k = 3 # Number of treatments

# Pre-process numerical covariates.
num_covariates = c("age_rand", "ht", "pwl_wt")

# Normalize numerical covariates to be between -1 and +1
Xn = trial.data[,num_covariates]
X.scale.n = apply(Xn, MARGIN = 2, FUN = function(x){a = (x-min(x))/(max(x)-min(x)); 2*a - 1})

# Pre-process categorical variables
cat_covariates = c("female")
Xc = trial.data[,cat_covariates]
Xc[Xc == 0] = -1

# Put covariate data together.
s.covariates = cbind("Intercept" = 1, X.scale.n, Xc)
H.full = as.matrix(s.covariates)
m = ncol(H.full) - 1 # Number of covariates.

# Create training and test sets.
set.seed(6836845)
id.select = sample(1:nrow(H.full), 115, replace = FALSE)
H = H.full[id.select,]
Zs = H.full[-id.select,]

# Save training matrix for genetic algorithm in Python
# write.csv(H, file = "Hmat.csv")

# Set actual clinical trial.
TRT.real = trial.data[id.select, "arm"]

## Read trial obtained from genetic algorithm.
X.gen.Iopt = read.csv("GA_asignacion_groups_Hmat_P100_G200_C1p0_M0p1_T2_E3.csv")[,-1]

# Evaluate genetic algorithm. 
subj.variances.genetic = matrix(NA, ncol = k, nrow = nrow(Zs))
for(j in 1:k) {
  x = (X.gen.Iopt == j)
  Hgroup = H[x,]
  R <- t(Hgroup)%*%Hgroup
  U <- solve(R)
  subj.variances.genetic[,j] = apply(Zs, MARGIN = 1, function(x) t(x)%*%U%*%x)
}

# Evaluate original trial.
X.original = rep(NA, length(TRT.real))
X.original[TRT.real == "1:LO"] = 1
X.original[TRT.real == "2:MOD"] = 2
X.original[TRT.real == "3:HI"] = 3
subj.variances.original = matrix(NA, ncol = k, nrow = nrow(Zs))
for(j in 1:k) {
  x = (X.original == j)
  Hgroup = H[x,]
  R <- t(Hgroup)%*%Hgroup
  U <- solve(R)
  subj.variances.original[,j] = apply(Zs, MARGIN = 1, function(x) t(x)%*%U%*%x)
}

zs.sort.original <- matrix(0, ncol = k, nrow = nrow(Zs))
zs.sort.gentic <- matrix(0, ncol = k, nrow = nrow(Zs))
for (i in 1:k){
  zs.sort.original[,i] <- sort(subj.variances.original[,i])
  zs.sort.gentic[,i] <- sort(subj.variances.genetic[,i])
}

colnames(zs.sort.original) = c("Group 1", "Group 2", "Group 3")
colnames(zs.sort.gentic) = c("Group 1", "Group 2", "Group 3")

mat_original <- data.frame(ID = 1:nrow(zs.sort.original), Method = "original", zs.sort.original)
mat_genetic  <- data.frame(ID = 1:nrow(zs.sort.gentic), Method = "genetic", zs.sort.gentic)

final_result <- rbind(mat_original, mat_genetic)


# Construct the FDS plot.
dat.plot = melt(final_result, c("Method", "ID"))
fds = ggplot(data = dat.plot, aes(x = ID/nrow(Zs), y = value, 
                                  color = variable, linetype = Method)) + geom_line(linewidth = 2.5)
fds = fds + xlab("Fraction of Subjects") + ylab("Relative GSP Variance") 
fds = fds + plot.theme + ylim(c(0, 0.4))
ggsave("figures/FDS.pdf", fds, width=11, height=8.5)

#final.data.Genetic = dat.plot
