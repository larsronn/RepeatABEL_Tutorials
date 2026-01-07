library(RepeatABEL)
rGLSinteraction <-
function(formula.FixedEffects = y ~ 1, genabel.data, phenotype.data, id.name="id", GRM = NULL, V = NULL, memory=1e8, interaction.term = NULL) {
	#Check input data
	if (class(genabel.data)!="gwaa.data2") stop("The input of genabel.data is not a GenABEL object")
	if (is.null(genabel.data@phdata$id)) stop("IDs not given as id in the phdata list")
    if (!is.null(GRM)) { if(!isSymmetric(GRM)) warning("The given GRM must be a symmetric matrix!") }
    if (!is.null(V)) { if(!isSymmetric(V)) warning("The given V must be a symmetric matrix!") }
    if (!is.null(interaction.term)) {
    	if (!(interaction.term%in%all.vars(formula.FixedEffects))) stop("interaction.term must be included in formula.FixedEffects")
    	print("If the interaction term is a factor, remember to define it as such before running this function.")
    }
    V.input  <- V
    #require(hglm)
    #require(GenABEL)
    #Get trait name
	trait <- all.vars(formula.FixedEffects)[1]
    #Remove NAs from phenotypic data
    y.all <- phenotype.data[,names(phenotype.data)%in%trait]
    phenotype.data <- phenotype.data[!is.na(y.all),]
	#Connect IDs in GenABEL data set with IDs in the phenotype file
	id1 <- phenotype.data[,names(phenotype.data)%in%id.name] #ID for phenotype data
	id2 <- genabel.data@phdata$id #ID for genotype data
    test1 <- id1%in%id2
    test2 <- id2%in%id1
    genabel.data <- keep_gwaa_data(genabel.data, which(test2)) #Exclude individuals having no phenotype information
    phenotype.data = phenotype.data[test1,] #Exclude individuals having no genotype information
    id1 <- phenotype.data[,names(phenotype.data)%in%id.name] #ID for phenotype data for cleaned data
	id2 <- genabel.data@phdata$id #ID for genotype data for cleaned data
	#####################
	#Construct incidence matrix for repeated observations
    N=length(id2)
    n=length(id1)
    indx <- numeric(n)
    for (i in 1:N) {
        indx <- indx+i*(id1%in%id2[i])
    }
    Z.indx <- diag(N)[indx,]
    #Construct response and design matrix
	y <- phenotype.data[,names(phenotype.data)%in%trait] #Create the response variable
    X <- model.matrix(formula.FixedEffects, data=phenotype.data) #Fixed effect design matrix
    if (is.null(interaction.term)) X.int <- 1
    if (!is.null(interaction.term)) X.int <-  model.matrix( as.formula(paste("~",interaction.term)), data=phenotype.data)
	#####################
 if (is.null(V)) {
    #Construct GRM
    if (is.null(GRM)) {
      autosomalMarkers <- which(chromosome(genabel.data)!= "X")
      GRM <- compute.GRM(genabel.data@gtdata@gtps[, snpnames(genabel.data)[autosomalMarkers]])
    }
    eig <- eigen(GRM)
    if (max(diag(GRM))>1.6) print("There seems to be highly inbred individuals in your data")
    if (min(eig$values < -0.5)) print("The genetic relationship matrix is far from positive definite")
    non_zero.eigenvalues = eig$values>(1e-6) #Put numerically small eigenvalues to zero
    eig$values[!non_zero.eigenvalues]=0
    print("GRM ready")
    #####################
    #Fit hglm
    Z.GRM <- ( eig$vectors%*%diag(sqrt(eig$values)) )[indx,]
    Z <- (cbind(Z.GRM, Z.indx))
    mod1 <- hglm(y=y, X=X, Z=Z, RandC = c(ncol(Z.GRM),ncol(Z.indx)), maxit = 200)
    if (mod1$Converge!="converged") stop("The variance component estimation did not converge in 200 iterations. Try to estimate them separately and provide the estimated (co)variance matrix V as input. \n\n")
    print("Variance component estimation ready")
    #####################
    #Construct rotation matrix
    ratio <- mod1$varRanef/mod1$varFix
    V <- constructV(Z=Z, RandC = c(ncol(Z.GRM),ncol(Z.indx)), ratio)
  }
	eig.V <- eigen(V)
	transf.matrix <- diag(1/sqrt(eig.V$values))%*%t(eig.V$vectors)
	y.new <- transf.matrix%*%y
	X.new <- transf.matrix%*%X
	print("Rotation matrix ready")
	#####################
	#Fit a linear model for each SNP
	SNP.matrix <- as.matrix(genabel.data@gtdata@gtps)
    if (sum(is.na(SNP.matrix)) > 0) {
        SNP.matrix <- SmoothSNPmatrix(SNP.matrix)
    }
	m = ncol(SNP.matrix)
	p.val <- SNP.est <- rep(1, m)
	colnames(X.new) <- as.character(1:ncol(X.new)) #To avoid columns having strange names
	print("Rotate LMM started")
    #Fit using QR factorization
    #Null model
	qr0 <- qr(X.new)
	est0 <- qr.coef(qr0,y.new)
	res <- y.new - X.new%*%est0
    n <-length(y.new)
	RSS.0 <- sum(res^2)/n
	#Split computations into reasonably sized blocks
    if (memory < n) memory <- n
	step.size <- floor(memory/n)
    steps <- ceiling(m/step.size)
	jj=1
	kk=0
	#########
	ColMult <- function(A,B) {
		mat <- matrix(0,nrow(A), ncol(A)*ncol(B))
		for (i in 1:ncol(A)){
			mat[,((i-1)*ncol(B)+1):(i*ncol(B))] <- A[,i]*B
		}
		return(mat)
	}
	#########
	for (step.i in 1:steps) {
		if (step.i==steps) kk=kk+m%%step.size else kk=kk+step.size
		markers.to.fit = jj:kk
		if (length(X.int)==1) snp.new <- transf.matrix%*%SNP.matrix[indx,markers.to.fit]
		if (length(X.int)>1) snp.new <- transf.matrix%*%ColMult(SNP.matrix[indx,markers.to.fit], X.int)
		mm=0
		for (j in markers.to.fit) {
			mm=mm+1
			if (length(X.int)==1) X1 <- cbind(snp.new[,mm], X.new)
			if (length(X.int)>1) X1 <- cbind(snp.new[,((mm-1)*ncol(X.int)+1):(mm*ncol(X.int))], X.new)
			qr1 <- qr(X1)
			est1 <- qr.coef(qr1,y.new)
			res <- y.new - X1%*%est1
			RSS.1 <- sum(res^2)/n
			SNP.est[j] <- est1[1]
			LRT <- -n*(log(RSS.1)-log(RSS.0))
			if (length(X.int)==1) {df <- 1} else {df <- ncol(X.int)}
			p.val[j] <-1 - pchisq(LRT, df=df)
        }
        jj=jj+step.size
	}
	
	print("Rotate LMM ready")
	#####################
	qt.results <- Create_gwaa_scan2(genabel.data, p.val, SNP.est)
	if (is.null(V.input)) qt.results@call$hglm <- mod1
	return(qt.results)
}

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

snpnames <- function(genabel.data) {
  return(genabel.data@gtdata@snpnames)
}
