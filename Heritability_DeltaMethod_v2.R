library(RepeatABEL)

# A help function to extract genotypes
keep_gwaa_data <- function(genabel.data, indx.keep = NULL) {
  #if (is.integer(indx.keep)) indx.keep = which(indx.keep)
  indx.keep <- as.numeric(indx.keep)
  genabel.data@gtdata@gtps = genabel.data@gtdata@gtps[indx.keep,]
  genabel.data@gtdata@nids = nrow(genabel.data@gtdata@gtps)
  genabel.data@gtdata@idnames = genabel.data@gtdata@idnames[indx.keep]
  tmp = as.data.frame(genabel.data@phdata[indx.keep,])
  names(tmp)=names(genabel.data@phdata)
  genabel.data@phdata = tmp
  return(genabel.data)
}

#Function to estimate h2 and its SE
get.SEh2 <- function(formula.FixedEffects = y ~ 1, genabel.data, phenotype.data, id.name = "id", GWAS.output, GRM = NULL) {
  if (class(genabel.data) != "gwaa.data2") 
    stop("The input of genabel.data is not a gwaa.data2 object")
  if (is.null(genabel.data@phdata$id)) 
    stop("IDs not given as id in the phdata list")
  if (is.null(GRM)) GRM <- compute.GRM(genabel.data@gtdata@gtps)
  trait <- all.vars(formula.FixedEffects)[1]
  y.all <- phenotype.data[, names(phenotype.data) %in% trait]
  phenotype.data <- phenotype.data[!is.na(y.all), ]
  id1 <- phenotype.data[, names(phenotype.data) %in% id.name]
  id2 <- genabel.data@phdata$id
  test1 <- id1 %in% id2
  test2 <- id2 %in% id1
  #genabel.data <- genabel.data[test2, ]
  genabel.data <- keep_gwaa_data(genabel.data, which(test2)) #Exclude individuals having no phenotype information
  phenotype.data <- phenotype.data[test1, ]
  id1 <- phenotype.data[, names(phenotype.data) %in% id.name]
  id2 <- genabel.data@phdata$id
  N = length(id2)
  n = length(id1)
  indx <- numeric(n)
  for (i in 1:N) {
    indx <- indx + i * (id1 %in% id2[i])
  }
  Z.indx <- diag(N)[indx, ]
  y <- phenotype.data[, names(phenotype.data) %in% trait]
  X <- model.matrix(formula.FixedEffects, data = phenotype.data)
  eig <- eigen(GRM)
  non_zero.eigenvalues <- eig$values > (1e-06)
  eig$values[!non_zero.eigenvalues] <- 0
  #print("GRM ready")
  Z.GRM <- (eig$vectors %*% diag(sqrt(eig$values)))[indx, ]
  Z <- cbind(Z.GRM, Z.indx)
  mod1b <- GWAS.output@call$hglm
  V <- constructV(Z=Z,RandC = c(ncol(Z.GRM), ncol(Z.indx)), ratio=mod1b$varRanef/mod1b$varFix)
  V <- mod1b$varFix*V #Bug fixed by LRN 24 May 2016
  tr <- function(X) sum(diag(X))
  invV <- solve(V)
  P <- invV - invV%*%X%*%solve( t(X)%*%invV%*%X )%*%t(X)%*%invV
  b <- c( mod1b$varRanef[2]+mod1b$varFix, -mod1b$varRanef[1], -mod1b$varRanef[1]) / ( (sum(mod1b$varRanef) + mod1b$varFix)^2 )
  row1 <- cbind( tr( P%*%tcrossprod(Z.GRM)%*%P%*%tcrossprod(Z.GRM)), tr( P%*%tcrossprod(Z.GRM)%*%P%*%tcrossprod(Z.indx)), tr( P%*%tcrossprod(Z.GRM)%*%P))
  row2 <- cbind( tr( P%*%tcrossprod(Z.indx)%*%P%*%tcrossprod(Z.GRM)), tr( P%*%tcrossprod(Z.indx)%*%P%*%tcrossprod(Z.indx)), tr( P%*%tcrossprod(Z.indx)%*%P))
  row3 <- cbind( tr( P%*%P%*%tcrossprod(Z.GRM)), tr( P%*%P%*%tcrossprod(Z.indx)), tr( P%*%P))
  C <- 0.5*rbind(row1, row2, row3)
  if (!isSymmetric(C)) warning("Non-symmetric matrix constructed")
  b <- matrix(b,1,3)
  #print(b)
  #print(round(cov2cor(solve(C)),3))
  h2 <- mod1b$varRanef[1]/(sum(mod1b$varRanef)+mod1b$varFix)
  return( c(h2, sqrt( b%*%solve(C)%*%t(b) ) ) ) 
}

###########################################################################################################################
################
#Example 1
set.seed(1234)
Gen.Data <- simulate_gendata(n=100, p=200)
Phen.Data <- simulate_PhenData(y ~ 1, genabel.data=Gen.Data,
                                n.obs=rep(4, nids(Gen.Data)), SNP.eff=2, SNP.nr=100, VC=c(1,1,1))
GWAS1 <- rGLS(y ~ 1, genabel.data = Gen.Data, phenotype.data = Phen.Data)
plot(GWAS1, main="")
summary(GWAS1)
#Summary for variance component estimation without SNP effects
summary(GWAS1@call$hglm)


h2.SE <- get.SEh2(formula.FixedEffects = y ~ 1, genabel.data=Gen.Data, phenotype.data=Phen.Data, GWAS.output=GWAS1) 
cat("The estimated heritability is ", round(h2.SE[1],3), "with an approximate SE of ", round(h2.SE[2],4), ".","\n" )
