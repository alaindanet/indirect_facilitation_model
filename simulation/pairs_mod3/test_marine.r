setwd("~/Bureau")
# data <- read.table("genotypes_f2.csv", header=T, sep=",",row.names=1)
data <- read.table("genotype_correction.txt", header=T, sep=",",row.names=1,)

freq <- matrix(data=NA, nrow=5, ncol=ncol(data), dimnames=list(c("E/E","E/G","G/G","E","G"),colnames(data)))
library(genetics)
for (i in 1:ncol(data)) {
	tempdata <- data[,i]; tempdata <- tempdata[which(tempdata!="NA/NA")]

	temp <- summary(genotype(tempdata))

	nameallele <- dimnames(temp$allele.freq)[[1]]
	if (length(nameallele)==2) {
		freq[c("G","E"),i] <- temp$allele.freq[c("G","E"),"Count"]
	} else {
		trickallele <- grep("G|E", dimnames(temp$genotype.freq)[[1]], perl=TRUE)
		freq[nameallele[trickallele],i] <- temp$allele.freq[nameallele[trickallele],"Count"]
	}
	print(i)
	namegeno <- dimnames(temp$genotype.freq)[[1]]
	trickHE <- grep("E/G|G/E", dimnames(temp$genotype.freq)[[1]], perl=TRUE)

	if (length(namegeno)==3) {
		freq[c("G/G","E/G","E/E"),i] <- temp$genotype.freq[c("G/G",namegeno[trickHE],"E/E"),"Count"]	
	} else {
		trickHO <- grep("G/G|E/E", dimnames(temp$genotype.freq)[[1]], perl=TRUE)
		freq[c(namegeno[trickHO],"E/G"),i] <- temp$genotype.freq[c(namegeno[trickHO],namegeno[trickHE]),"Count"]
	}
}

# Pour le tableau de contingence, remplacer les NA par des 0
freq[which(is.na(freq))] <- 0
test <- matrix(c(freq[1:3,1],0.25,0.5,0.25),ncol=3, nrow=2, dimnames=list(c("obs","th"), c("E/E","E/G","G/G")), byrow=T)
fisher.test(test)

# Fisher's Exact Test for Count Data
contingence <- matrix(NA,ncol=3, nrow=2, dimnames=list(c("obs","th"), c("E/E","E/G","G/G")), byrow=T)

resfisher <- matrix(NA, ncol=ncol(freq), nrow=1, dimnames=list(paste("Pval",round(1-0.05/ncol(freq),3),sep=""), colnames(freq)))

for (i in 1:ncol(freq)) {
	contingence["obs",c("E/E","E/G","G/G")] <- freq[c("E/E","E/G","G/G"),i]
	ntot <- sum(contingence["obs",])
	contingence["th",c("E/E","E/G","G/G")] <- round(c(0.25*ntot,0.50*ntot, 0.25*ntot))
	resfisher[,i] <- fisher.test(contingence, conf.level=0.95,3), hybrid=T)$p.value 
}
resfisher[,which(resfisher<0.05)]# Classic significant threhold
resfisher[,which(resfisher<0.05/ncol(resfisher))] # Badass Boneferroni correction
write.table(resfisher,file="fisher.result.csv",sep=",")
