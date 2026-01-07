library(RepeatABEL)
source("rGLSinteraction.R")
set.seed(12345)
gen.data <- simulate_gendata(n=500, p=2000)
VC.poly <- VC.perm <- 0
VC.res <- 1
n.obs <- rep(4, nids(gen.data))
Phen.Sim <- simulate_PhenData(y ~ 1, 
genabel.data=gen.data, n.obs=n.obs, SNP.eff=0, 
SNP.nr=1000, VC=c(VC.poly,VC.perm,VC.res))
id.name="id"
phenotype.data<- Phen.Sim
genabel.data <- gen.data
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
    
SNP <- gen.data@gtdata@gtps[indx, 1000] 
table(SNP)
tmp <- sample(c(-1,1), 100, replace=TRUE)
####### PRODUCE YEAR EFFECTS FOR EACH INDIVIDUAL
#sd.year=10 #Standard Deviation of Year Effects
#beta <- rnorm(max(n.obs), 0, sd.year) #Simulated Year Effects
#beta <- 0.45*c(1,-1,1,-1) #A moderate effect
beta <- 2*c(1,-1,1,-1) 
year.effects <- years <- NULL
for (i in 1:sum(n.obs)) {
  yr.i <- i%%4
  if (yr.i==0) yr.i=4
  years <- c(years, yr.i)
	year.effects <- c(year.effects, beta[yr.i]*SNP[i])

}
########################
#### A FUNCTION TO ADD A VARIABLE TO A LIST
add.var <- function(x, add.new, new.name) {
  x[[length(x)+1]] <- add.new
	names(x)[length(x)] <- new.name
	return(x)
}
#######################
#ADDS A NEW PHENOTYPE WITH YEAR EFFECTS ADDED
Phen.Sim <- add.var(Phen.Sim, Phen.Sim$y + year.effects, "y.yrs")  
#ADD YEAR AS FACTOR
Phen.Sim <- add.var(Phen.Sim, as.factor(years), "YearsAll") 
#######################
Phen.Sim$Years <- as.numeric(Phen.Sim$Years )
lm0 <- lm(y.yrs~SNP+YearsAll, data=Phen.Sim)
summary(lm0)
lm1 <- lm(y.yrs~SNP+SNP*YearsAll, data=Phen.Sim)
summary(lm1)


#RUN THE ANALYSIS
GWAS.sim1 <- rGLS(y.yrs ~ YearsAll, genabel.data = gen.data, 
                  phenotype.data = Phen.Sim)
plot(GWAS.sim1)

GWAS.sim2 <- rGLSinteraction(y.yrs ~ YearsAll, genabel.data = gen.data, 
                             phenotype.data = Phen.Sim, interaction.term="YearsAll")
plot(GWAS.sim2)

#Genomic deflation in these simulated data??
estlambda(GWAS.sim1@results$P1df, method="median")
estlambda(GWAS.sim2@results$P1df, method="median")

#### What does the SNP effect look like for the different years
library(ggplot2)
ggplot(Phen.Sim, aes(x=factor(SNP), y=y.yrs)) +
  geom_boxplot() +
  facet_wrap(~YearsAll) +
  ggtitle("SNP effect per year")
