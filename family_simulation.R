# Generate a child from two parents (father=number of column in X (VCF without header) for the father, mother=number of column in VCF for the mother)

offspring <- function(X, father=NULL, mother=NULL) {
  fatherGT <- X[,father]
  motherGT <- X[,mother]
  fatherAll <- str_split_fixed(fatherGT, "/", 2)
  motherAll <- str_split_fixed(motherGT, "/", 2)
  snps <- nrow(X)
  f <- vector(length = snps)
  m <- vector(length = snps)
  f[1] <- sample(1:2, prob = c(0.5,0.5))[1]
  m[1] <- sample(1:2, prob = c(0.5,0.5))[1]
	for(n in 2:snps) {
	  d <- abs((X[n,2]-X[n-1,2])/1000)
	  if(d==0) {
	    prec <- 0
	  }else{
	    prec <- 0.5*(1-exp(-2*d))
	  }
	  if(f[n-1]==1) {
	    f[n] <- sample(1:2, prob = c(1-prec,prec))[1]
	  }else{
	    f[n] <- sample(1:2, prob = c(prec,1-prec))[1]
	  }
	  if(m[n-1]==1) {
	    m[n] <- sample(1:2, prob = c(1-prec,prec))[1]
	  } else {
	    m[n] <- sample(1:2, prob = c(prec,1-prec))[1]
		}
	}
  
  child <- sapply(1:snps, function(i) {
    paste(as.matrix(fatherAll[i, f[i]]), as.matrix(motherAll[i, m[i]]), 
          sep = "/")
  })
  
  X <- cbind(X, child)
  colnames(X) <- 1:ncol(X)
  return(X)

}


# Import the data from the founders

library(data.table)
library(stringr)

founders <- data.frame(fread("founders_sort.txt"))
colnames(founders) <- 1:ncol(founders)

# Run only one pedigree type

#type 1 pedigree (16 ind)

family <- offspring(founders, 10, 11)
family <- offspring(family, 12, 13)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)
family <- offspring(family, 14, 15)

#type 2 pedigree (16 ind)

family <- offspring(founders, 10, 11)
family <- offspring(family, 10, 11)
family <- offspring(family, 10, 11)
family <- offspring(family, 10, 11)
family <- offspring(family, 12, 14)
family <- offspring(family, 12, 14)
family <- offspring(family, 12, 14)
family <- offspring(family, 12, 14)
family <- offspring(family, 13, 14)
family <- offspring(family, 13, 14)
family <- offspring(family, 13, 14)
family <- offspring(family, 13, 14)

#type 3 pedigree (16 ind)

family <- offspring(founders, 10, 11)
family <- offspring(family, 10, 11)
family <- offspring(family, 14, 12)
family <- offspring(family, 14, 12)
family <- offspring(family, 14, 12)
family <- offspring(family, 14, 12)
family <- offspring(family, 14, 12)
family <- offspring(family, 13, 15)
family <- offspring(family, 13, 15)
family <- offspring(family, 13, 15)
family <- offspring(family, 13, 15)
family <- offspring(family, 13, 15)

#type 4 pedigree (12 ind)

family <- offspring(founders, 10, 11)
family <- offspring(family, 11, 12)
family <- offspring(family, 11, 12)
family <- offspring(family, 13, 17)
family <- offspring(family, 14, 18)
family <- offspring(family, 15, 19)


# Save simulated family

write.table(family, file = "simulated_family.txt", quote = F, sep = "\t", row.names = F, col.names = F)














