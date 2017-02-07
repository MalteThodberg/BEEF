library(polyester)
library(ABCdata)
library(edgeR)

### SEED

Seed <- 2015
set.seed(Seed)

### Paramter

# Simulate based on rat neuron data
data(RatNeurons)

# Filter out low counts
filter <- apply(RatNeurons$Expression, 1, function(x) length(x[x>0])>=2)
counts <- RatNeurons$Expression[filter,]

# Get parameter for later
params <-  get_params(counts)

### Design

design <- data.frame(bio=c(replicate(50, "x"),
													 replicate(50, "y")))
design$batch <- LETTERS[1:4]

mod <- model.matrix(~0+batch+bio, data=design)

modI <- model.matrix(~batch+bio, data=design)
mod0 <- modI[,c(1, ncol(modI))]
modB <- modI[,2:4]
### Effect sizes

n_genes <- 10000

effect_sizes <- list(batchA=0.10,
										 batchB=0.10,
										 batchC=0.10,
										 batchD=0.10,
										 bioy=0.05)

betas <- sapply(effect_sizes, function(x) sample(c(rnorm(n_genes*x), rep(0, n_genes-(n_genes*x)))))

### Actual simulation
sim_dat <- create_read_numbers(params$mu,
															 params$fit,
															 params$p0,
															 #m=dim(counts)[1],
															 #n=dim(counts)[2],
															 beta=betas,
															 mod=mod,
															 seed=Seed)

# Add colum and row names

### Normalize

# Normalizse
dge <- DGEList(subset(sim_dat, rowSums(cpm(sim_dat) >= 1) >= 3))
dge <- calcNormFactors(dge, method="upperquartile")
em <- cpm(dge, log=TRUE, pseudo.count=1)

# Rename
colnames(em) <- paste0("Sample", 1:ncol(em))
rownames(em) <- paste0("Gene", 1:nrow(em))

### Clean and save
batch4class2 <- list(x=t(sim_dat), y=factor(design$bio))

devtools::use_data(batch4class2, overwrite=TRUE)
