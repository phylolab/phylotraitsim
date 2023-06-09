library(stringr)
library(ape)

args <- commandArgs(trailingOnly=TRUE)
directory <- args[1]
#directory <- "."
seed <- args[2]
#seed <- 123
specType <- 4
landscape <- 1

data <- read.csv(paste0(directory,"/test_numInd_perArea_perSpecies_2_",seed,".csv"),h=F)
treeAll <- read.tree(file=paste0(directory,"/output_newickTreeAll_a_",specType,"_f_",landscape,"_",seed,".tre"))
landscape_table <- read.csv(paste0(directory,"/landscapes_per_area_2.csv"),header=TRUE)
K_areas <- landscape_table$Ka

lenData <- length(data[,1])
#numSpecies <- 103
numSpecies <- length(treeAll$tip.label)
numAreas <- length(K_areas)
dred <- data[3:lenData,1:(numSpecies*numAreas)]
#dred[dred>0] <- 1

print("Data read.")

inds_areas <- matrix(0,nrow=length(dred[,1]),ncol=numAreas)
species_areas <- numeric(numAreas)
species_areas_2 <- numeric(numAreas)
species_areas_3 <- numeric(numAreas)
species_areas_4 <- numeric(numAreas)
currentSpecies <- dred[(lenData-2),]
currentSpecies_2 <- dred[30000,]
currentSpecies_3 <- dred[15000,]
currentSpecies_4 <- dred[7500,]
currentSpecies[currentSpecies>0] <- 1
currentSpecies_2[currentSpecies_2>0] <- 1
currentSpecies_3[currentSpecies_3>0] <- 1
currentSpecies_4[currentSpecies_4>0] <- 1
for (i in 1:numAreas) {
    inds_areas[,i] <- apply(dred[,seq(i,numSpecies*numAreas,numAreas)],1,sum)
    species_areas[i] <- sum(currentSpecies[seq(i,numSpecies*numAreas,numAreas)])
    species_areas_2[i] <- sum(currentSpecies_2[seq(i,numSpecies*numAreas,numAreas)])
    species_areas_3[i] <- sum(currentSpecies_3[seq(i,numSpecies*numAreas,numAreas)])
    species_areas_4[i] <- sum(currentSpecies_4[seq(i,numSpecies*numAreas,numAreas)])
}

print("Done!")

write.table(inds_areas,file=paste0(directory,"/res_analyze_K_",seed,".csv"),sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(t(c(species_areas,species_areas_2,species_areas_3,species_areas_4)),file=paste0(directory,"/res_species_",seed,".csv"),sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)


#par(mfrow=c(1,2))
#plot(1:length(inds_areas[,1]),inds_areas[,1],type="l",xlab="Generations",ylab="N",ylim=c(0,5000))
#for (i in 2:numAreas) {
#    lines(1:length(inds_areas[,1]),inds_areas[,i],col=i)
#}
#legend("topleft",legend=c(1:numAreas),pch=15,col=1:numAreas,bty="n")

#barplot(apply(inds_areas,2,sum),names=1:numAreas,ylab="N",xlab="Areas")
