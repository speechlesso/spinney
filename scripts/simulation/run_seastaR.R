args <- commandArgs(trailingOnly = TRUE)

if(length(args) != 3) {
  stop("Please provide two inputs!")
}

sptree_path <- args[1]
outfile <- args[2]
sigma2 <- as.numeric(args[3])

# load library
.libPaths('~/Rlibs')
library(phytools)
library(seastaR)
sptree <- seastaR::parse_input_file(sptree_path, genetrees = FALSE)

Cstar_matrix <- seastaR::get_full_matrix(sptree)

test_trait <- seastaR::simulate_traits(100,Cstar_matrix, sigma2)

# write output
write.csv(test_trait, outfile, quote = FALSE, row.names=FALSE)

# inference

inferred_rate = as.matrix(0L, nrow = dim(test_trait)[1])
for (nrow in 1:dim(test_trait)[1])
  {
    inferred_rate[nrow] = seastaR::sigma2_inference(Cstar_matrix, test_trait[nrow,])
  }

print(c(sigma2, mean(inferred_rate)))
print(sd(inferred_rate))
  
