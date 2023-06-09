library(stringr)
library(ape)

args <- commandArgs(trailingOnly=TRUE)
directory <- args[1]
seed <- args[2]
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
dred[dred>0] <- 1

print("Data read.")

speciationsPerArea <- numeric(numAreas)
for (j in seq(1,(numSpecies*numAreas),numAreas)) { #Loop per species
    for (i in 1:lenData) { #Loop per generation
        if (sum(dred[i,j:(j+numAreas-1)])!=0) {
            speciationsPerArea[which(dred[i,j:(j+numAreas-1)]==1)] <- speciationsPerArea[which(dred[i,j:(j+numAreas-1)]==1)]+1
            break
        }
    }
}

print("Speciations done.")

extinctionsPerArea <- numeric(numAreas)
extinctionRatePerArea <- numeric(numAreas)
speciationRatePerArea <- numeric(numAreas)
totSum <- numeric(numAreas)
for (area in 1:numAreas) {
    for (a in seq(area,(numSpecies*numAreas),numAreas)) {
        ha <- paste(dred[,a],collapse="")
        extinctionsPerAreaPerSpecies <- dim(str_locate_all(ha,"10")[[1]])[1]
        extinctionsPerArea[area] <- extinctionsPerArea[area] + extinctionsPerAreaPerSpecies
        #extinctionRatePerArea[area] <- extinctionRatePerArea[area] + extinctionsPerAreaPerSpecies/sum(dred[,a])
        totSum[area] <- totSum[area] + sum(dred[,a]) #Sums the total time spent in area 'area'.
    }
    extinctionRatePerArea[area] <- extinctionsPerArea[area]/totSum[area]
    speciationRatePerArea[area] <- speciationsPerArea[area]/totSum[area]
}

print("Extinctions done.")
write.table(t(c(speciationsPerArea,extinctionsPerArea,totSum)),file=paste0(directory,"/res_analyze_full_",seed,".csv"),sep=",",quote=FALSE,row.names=FALSE,col.names=FALSE)

#par(mfrow=c(2,2))
#barplot(speciationsPerArea,xlab="Area",ylab="Number of speciations",names.arg=1:numAreas)
#barplot(extinctionsPerArea,xlab="Area",ylab="Number of (local) extinctions",names.arg=1:numAreas)
#barplot(speciationRatePerArea,xlab="Area",ylab="Speciation rate",names.arg=1:numAreas)
#barplot(extinctionRatePerArea,xlab="Area",ylab="Extinction rate",names.arg=1:numAreas)
