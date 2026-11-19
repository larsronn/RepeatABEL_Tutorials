

rGLS_dom <- function (formula.FixedEffects = y ~ 1, genabel.data, phenotype.data,
                      id.name = "id", GRM = NULL, V = NULL, memory = 1e+08, verbose = TRUE)
{
  if (!is(genabel.data, "gwaa.data2"))
    stop("The input of genabel.data is not a gwaa.data2 object")
  if (is.null(genabel.data@phdata$id))
    stop("IDs not given as id in the phdata list")
  if (!is.null(GRM)) {
    if (!isSymmetric(GRM))
      warning("The given GRM must be a symmetric matrix!")
  }
  if (!is.null(V)) {
    if (!isSymmetric(V))
      warning("The given V must be a symmetric matrix!")
  }
  V.input <- V
  trait <- all.vars(formula.FixedEffects)[1]
  y.all <- phenotype.data[, names(phenotype.data) %in% trait]
  phenotype.data <- phenotype.data[!is.na(y.all), ]
  id1 <- phenotype.data[, names(phenotype.data) %in% id.name]
  id2 <- genabel.data@phdata$id
  test1 <- id1 %in% id2
  test2 <- id2 %in% id1
  genabel.data <- RepeatABEL:::keep_gwaa_data(genabel.data, which(test2))
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
  if (is.null(V)) {
    if (is.null(GRM)) {
      autosomalMarkers <- which(chromosome(genabel.data) !=
                                  "X")
      GRM <- compute.GRM(genabel.data@gtdata@gtps[, RepeatABEL:::snpnames(genabel.data)[autosomalMarkers]])
    }
    eig <- eigen(GRM)
    if (max(diag(GRM)) > 1.6)
      warning("There seems to be highly inbred individuals in your data")
    if (min(eig$values < -0.5))
      warning("The genetic relationship matrix is far from positive definite")
    non_zero.eigenvalues <- eig$values > (1e-06)
    eig$values[!non_zero.eigenvalues] <- 0
    if (verbose)
      cat("GRM ready \n")
    Z.GRM <- (eig$vectors %*% diag(sqrt(eig$values)))[indx,
    ]
    Z <- (cbind(Z.GRM, Z.indx))
    mod1 <- hglm(y = y, X = X, Z = Z, RandC = c(ncol(Z.GRM),
                                                ncol(Z.indx)), maxit = 200)
    if (mod1$Converge != "converged")
      stop("The variance component estimation did not converge in 200 iterations. Try to estimate them separately and provide the estimated (co)variance matrix V as input. \n\n")
    if (verbose)
      cat("Variance component estimation ready \n")
    ratio <- mod1$varRanef/mod1$varFix
    V <- constructV(Z = Z, RandC = c(ncol(Z.GRM), ncol(Z.indx)),
                    ratio)
  }
  eig.V <- eigen(V)
  transf.matrix <- diag(1/sqrt(eig.V$values)) %*% t(eig.V$vectors)
  y.new <- transf.matrix %*% y
  X.new <- transf.matrix %*% X
  if (verbose)
    cat("Rotation matrix ready \n")
  SNP.matrix <- as.matrix(genabel.data@gtdata@gtps)
  if (sum(is.na(SNP.matrix)) > 0) {
    SNP.matrix <- RepeatABEL:::SmoothSNPmatrix(SNP.matrix)
  }
  m <- ncol(SNP.matrix)
  p.val <- SNP.est <- SNP.est.d <- rep(1, m)
  colnames(X.new) <- as.character(1:ncol(X.new))
  if (verbose)
    cat("Rotate LMM started \n")
  qr0 <- qr(X.new)
  est0 <- qr.coef(qr0, y.new)
  res <- y.new - X.new %*% est0
  n <- length(y.new)
  RSS.0 <- sum(res^2)/n
  if (memory < n)
    memory <- n
  step.size <- floor(memory/n)
  steps <- ceiling(m/step.size)
  jj = 1
  kk = 0
  for (step.i in 1:steps) {
    if (step.i == steps)
      kk <- kk + m%%step.size
    else kk <- kk + step.size
    markers.to.fit <- jj:kk
    snp.new <- transf.matrix %*% SNP.matrix[indx, markers.to.fit]
    ## Heterozygote coding for dominance coefficient
    snp.het.coded <- matrix(as.numeric(SNP.matrix[indx, markers.to.fit] == 1),
                            ncol = ncol(snp.new), nrow = nrow(snp.new))
    snp.d.new <- transf.matrix %*% snp.het.coded
    mm <- 0
    for (j in markers.to.fit) {
      mm <- mm + 1
      X1 <- cbind(snp.new[, mm], snp.d.new[, mm], X.new)
      qr1 <- qr(X1)
      est1 <- qr.coef(qr1, y.new)
      res <- y.new - X1 %*% est1
      RSS.1 <- sum(res^2)/n
      SNP.est[j] <- est1[1]
      ## Add estimates for dominance coefficient
      SNP.est.d[j] <- est1[2]
      LRT <- -n * (log(RSS.1) - log(RSS.0))
      ## Changed degrees of freedom
      p.val[j] <- 1 - pchisq(LRT, df = 2)
    }
    jj <- jj + step.size
  }
  if (verbose)
    cat("Rotate LMM ready \n")
  ## Modify results data structure to include dominance estimates
  results_df <- data.frame(SNP = genabel.data@gtdata@snpnames,
                           estA = SNP.est, estD = SNP.est.d, p.val = p.val)
  if (is.null(V.input)) {
    results <- list(results = results_df, hglm = mod1)
  } else {
    results <- results_df
  }

  return(results)
}

