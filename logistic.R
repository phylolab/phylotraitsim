args <- commandArgs(trailingOnly=TRUE)
seed <- args[1]

logistic <- function(x,x_mid,k,x_max){
  return(x_max/(1+exp(-k*(x-x_mid))))
}
gen_midpoint_trajectory <- function(numGenerations,evo.rate) {
  traj <- numeric(numGenerations)
  min_midPoint <- 2 #Area 1+1
  max_midPoint <- 8 #Area 10-1
  traj[1]=5 #Stat in the middle area
  for (i in 2:numGenerations) {
    x=rnorm(1,mean=traj[i-1],sd=evo.rate)
    if (x>max_midPoint) { #Checks for boundaries
      x=max_midPoint-(x-max_midPoint)
    }
    if (x<min_midPoint) { #Checks for boundaries
      x=min_midPoint+(abs(x-min_midPoint))
    }
    traj[i]=x
  }
  return(traj)
}
numGenerations <- 25000
numAreas <- 9
maxOptimum <- 5 #x_max
minOptimum <- 0
steepness <- 2 #k
evo.rate <- as.numeric(args[2])
#midpoint.traj <- rep(5,numGenerations) #Area number where the midpoint of the logistic is located.
#midpoint.traj <- seq(2,8,length.out=numGenerations)
midpoint.traj <- gen_midpoint_trajectory(numGenerations,evo.rate)
write.table(as.data.frame(midpoint.traj),file=paste0("midpoint_trajectory_",seed,".csv"),append=FALSE,quote=FALSE,sep=",",row.names = FALSE,col.names = FALSE)
#par(mfrow=c(1,2))
#plot(midpoint.traj,type="l",ylim=c(1,numAreas))

for (j in 1:length(midpoint.traj)) {
  midpoint <- midpoint.traj[j]
  A <- numeric(numAreas)
  for (i in 1:numAreas) {
    A[i] <- logistic(i,midpoint,steepness,maxOptimum)
  }
#  if (j==1) {
#    plot(A,xlab="Areas",ylab="Optimum")
#  } else {
#    points(A)
#  }
  if (j==1) {
    write.table(t(as.data.frame(round(A,3))),file=paste0("optima_trajectories_movingMidpoint_",seed,".csv"),append=FALSE,quote=FALSE,sep=",",row.names = FALSE,col.names = FALSE)
  } else {
    write.table(t(as.data.frame(round(A,3))),file=paste0("optima_trajectories_movingMidpoint_",seed,".csv"),append=TRUE,quote=FALSE,sep=",",row.names = FALSE,col.names = FALSE)
  }
}
