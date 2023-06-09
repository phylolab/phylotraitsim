###################-FORWARD SIMULATOR OF PHENOTYPES-##########################
###################-Written by Pablo Duchen and Daniele Silvestro-############

#set.seed(1234)

library(ape)
library(geiger)
library(optparse)
library(Rcpp)
library(data.table)

#Test 1 vs 4 vs 8 vs unlimited (maybe use stepping stone with 100 islands) areas. 
#Run 1: Test different means, sd's and Ks per area.
#Initialize simulation with a species that has a high migration rate.
#Rampal Etienne papers on speciation. Add competition. Moving optima with time (Hansen). Study changes in extinction based on constant optima versus moving optima. Make migration matrix change with time. Also, implement specialists vs generalists (the landscape should be the product of the species landscape and the area landscape). Line of least resistance (Dolph). Or make sigma2_G more variable and inheritable to newer species (this for specialist vs generalist).

useCPP <- 0
maxNumExtantSpeciesAllowed <- 1000
maxNumSpeciesAllowed <- 250
plot_every_x_gen <- 1000
speciationPerArea <- 1 #If 0 speciation happens in all areas, regardles of geography.
initialArea <- sample(1:9,1) #CHANGE HERE THE NUMBER OF AREAS
Prob_sympatric <- function() return(0.5) #Prob allopatric is 1 minus Prob_sympatric.
get_sigma2_s <- function() return(runif(1,0.5,5)*10000) #This is the intra-specific sd2 optimum sigma2_s.
get_growth_rate <- function() return(runif(1,1.1,2)) #This is the growth rate g. Alternatively use rbeta(100,1.2,5)+1 to center the draws around 1.2
get_sigma2_G <- function() return(runif(1,0.005,0.02)) #This is the parent-offspring variance Sigma2_G.
#get_P_migration_species_specific <- function() return(min(rexp(1,1/0.0001),0.999))
get_P_migration_species_specific <- function() return(sample(c(0.00001,0.0001,0.01),1,prob=c(0.10,0.70,0.20)))
get_species_phen_range <- function() return(runif(1,0,5))
get_K <- function(lenNewSpecies) { #This is K per species.
    max_species_K <- 4000
    if (lenNewSpecies>=max_species_K) {
        return(lenNewSpecies)
    } else {
        return(runif(1,lenNewSpecies,max_species_K))
    }
}
convert_Ks <- function(K_areas) { #The final K per species will also depend on the size of the area.
    N <- sum(K_areas)
    ns <- K_areas/N
    f <- 1/min(ns)
    return(ns*f)
}

allColors <- rainbow(100)
colors <- allColors[1]

######-FUNCTIONS-######

fitness_density <- function(x,a=0.99,K=4,d=0){ #Function to choose between fitness landscapes
    if (landscape==0) {
        lik <- rep(0.01,length(x)) #Flat
    } else if (landscape==1) {
        lik <- dnorm(x,mean=trait_optimum,sd=sd_optimum) #One peak
    } else if (landscape==2) {
        depth <- 10; f <- 10; c <- 0.1;
        lik <- (depth*(x-0.5*max_trait_val)^2+f)*exp(-c*(x-0.5*max_trait_val)^2)-(max_trait_val*2) #Two peaks
    } else if (landscape==3) {
        lik = log( sqrt((1+a^2) - 2*a*cos( ((x+d*K)) * 1/K*(2*pi))  )) #Multiple peaks
    }  
    return (lik)
}

get_joint_mu_s2 <- function (mu_f,s2_f,mu_g,s2_g){ #Input vector of moms and vector of dads
    s2_fg = (s2_f*s2_g)/(s2_f+s2_g)
    mu_fg = (mu_f/s2_f + mu_g/s2_g) * s2_fg
    return(list(mu_fg,s2_fg))
    #return(c(mu_fg))
}

get_species_trait_values_unsorted <- function(current_pop,current_spec_ind){
    species_trait_val = c()
    species_trait_var = c()
    species_index =     c()
    count_species_ind = c()
    species_index_in_vector = NULL
    
    COUNTER = 1
    for (sp in unique(current_spec_ind)){
        m = current_pop[current_spec_ind==sp]
        species_trait_val =c(species_trait_val ,mean(m))
        species_trait_var =  c(species_trait_var ,var(current_pop[current_spec_ind==sp]))
        species_index = c( species_index,sp)
        count_species_ind = c(count_species_ind,length(m))
        species_index_in_vector[current_spec_ind==sp] = COUNTER
        COUNTER = COUNTER +1 
    }
    return( list(species_trait_val,species_trait_var,species_index,count_species_ind,species_index_in_vector) )
}

gen_couples <- function(current_pop_j,current_geography_j) {
    genSexes <- sample(c(1,2),length(current_pop_j),replace=TRUE) #We randomly assign sexes male or female to the individuals in species j.
    moms <- current_pop_j[which(genSexes==1)]
    moms_positions <- current_geography_j[which(genSexes==1)]
    M <- data.frame(moms,moms_positions)
    #M <- M[order(M$moms_positions,M$moms),] #order by geographic position and then by phenotype.
    M <- M[order(M$moms_positions),] #order by geographic position only.
    dads <- current_pop_j[which(genSexes==2)]
    dads_positions <- current_geography_j[which(genSexes==2)]
    D <- data.frame(dads,dads_positions)
    #D <- D[order(D$dads_positions,D$dads),] #order by geographic position and then by phenotype.
    D <- D[order(D$dads_positions),] #order by geographic position only.

    moms_sorted <- c()
    dads_sorted <- c()
    couples_positions <- c() #These are actually the couples locations (which areas are they in)
    for (i in unique(sort(current_geography_j))) { #This loop makes sure the couples are formed only within an area.
        M_area <- subset(M,M[,2]==i)
        D_area <- subset(D,D[,2]==i)
        lM <- length(M_area[,2])
        lD <- length(D_area[,2])
        m <- min(lD,lM)
        if (m<10) { #make sure there are at least 10 couples per area.
            next
        }
        moms_sorted <- c(moms_sorted,M_area[1:m,1])
        dads_sorted <- c(dads_sorted,D_area[1:m,1])
        couples_positions <- c(couples_positions,rep(i,m))
    }
    minl <- length(moms_sorted)

    return(list(moms_sorted,dads_sorted,minl,couples_positions))
}

sample_geo_individual <- function(i, migration_matrix,P_migration_species_specific){ #i is a vector of size 2 representing an area and the species id of a given individual.
    if (runif(1)<P_migration_species_specific[i[2]]) {
        areas <- 1:dim(migration_matrix)[1]
        return( sample(areas, 1, p=migration_matrix[i[1],]) ) #Make more efficient by deciding a priori which individuals will move, and use the apply only on those.
    } else {
        return(i[1])
    }
}

select_island <- function(x,islands,probMigration) { #Here x has to be a single element.
    w <- runif(1)
    if (w>probMigration) {
        return(x)
    } else {
        return(sample(islands,1))
    }
}

geography <- function(new_geography,new_pop_spec_ind,probMigration,migration_matrix,P_migration_species_specific,typeGeography) {
    #islands <- seq(-ceiling(numIslands/2),ceiling(numIslands/2))
    if (typeGeography==1) { #stepping stone
        migrants <- sample(c(0,1,-1),size=length(new_geography),replace=TRUE,prob=c((1-2*probMigration),probMigration,probMigration))
        new_geography <- new_geography+migrants
        new_geography[new_geography>max(islands)] <- max(islands) #These steps make sure the limit regions are not surpassed.
        new_geography[new_geography<min(islands)] <- min(islands)
    } else if (typeGeography==2) { #island model
        new_geography <- sapply(new_geography,select_island,islands,probMigration) #This line replaces each element of new_geography with a new region or island with probability probMigration
    } else if (typeGeography==3) { #General migration matrix
        new_geography <- apply(rbind(new_geography,new_pop_spec_ind), 2, FUN=sample_geo_individual, migration_matrix, P_migration_species_specific) 
    }
    return(new_geography)
}


migrationMat <- function(numIslands,probMigration,typeGeography,migration_table) {
    migration_matrix <- matrix(rep(0,numIslands^2),nrow=numIslands,byrow=T)
    if (typeGeography==2) { #island model, but where it is more likely to stay where you are.
        probs2 <- (probMigration)/(dim(migration_matrix)[1]-1) #the probability of migration is distributed among other areas.
        migration_matrix[,] <- probs2
        diag(migration_matrix) <- 1-probMigration
    } else if (typeGeography==1) { #stepping stone
        probs2 <- (probMigration)/2 #the probability of migration is distributed only among neighboring areas.
        diag(migration_matrix) <- 1-probMigration
        for (i in 1:numIslands) {
            if (i==1) {
                migration_matrix[i,i+1] <- probMigration
            } else if (i==numIslands) {
                migration_matrix[i,i-1] <- probMigration
            } else {
                migration_matrix[i,i+1] <- probs2
                migration_matrix[i,i-1] <- probs2
            }
        }
    } else if (typeGeography==3) { #Migration matrix provided by the user.
        migration_matrix <- migration_table
        migration_matrix <- as.matrix(migration_matrix)
    }
    return(migration_matrix)
}

assignParameters2species <- function(new_parameters_vec,current_spec_ind) { #alternative: new_parameters_vec[current_spec_ind]
    len <- length(current_spec_ind)
    res <- numeric(len)
    for (i in 1:len) {
        res[i] <- new_parameters_vec[current_spec_ind[i]]
    }
    return(res)
}

assign_death_prob_K <- function(species_id,area_id,growth_rate,K,free_space_area) {
    #n <- length(subset(new_geo_pop_spec,new_geo_pop_spec$new_geography==area_id & new_geo_pop_spec$new_pop_spec_ind==species_id)[,1])
    return( (growth_rate[species_id]-1)*n[species_id,area_id]/min(K[species_id],free_space_area[area_id]) )
}

Rcpp::cppFunction("
     std::vector<double> myDiag(Rcpp::NumericMatrix M) {  
     int ncolumns = M.ncol();
     std::vector<double> d(ncolumns);
     for (int i = 0; i < ncolumns; i++) {
             d[i] = M(i,i);
     }
     return d;
    }")

Rcpp::cppFunction("
     std::vector<double> assign_death_prob_K_cpp(Rcpp::NumericMatrix M_K,Rcpp::NumericMatrix M_IA) {  
     int ncolumns = M_IA.ncol();
     std::vector<double> dk(ncolumns);
     for (int i = 0; i < ncolumns; i++) {
         dk[i] = M_K((M_IA(1,i)-1),(M_IA(0,i)-1));
     }
     return dk;
    }")


writeEffectiveMigrations <- function(i,from_vec,to_vec,spec_ind) {
    #survived_geography_startingPoint <- backup_new_geography[survived_new_pop_ind] #Compare this vector with survived_geography and find out which ones are different and where they originated. Then match it with the species vector.
    indicesMoving <- which(from_vec!=to_vec)
    to_output <- c(length(indicesMoving),from_vec[indicesMoving],to_vec[indicesMoving],spec_ind[indicesMoving])
    
    write("\\",outfile_migration,append=TRUE)
    write(paste("Generation",i),outfile_migration,append=TRUE)
    write(to_output,outfile_migration,ncolumns=length(to_output),append=TRUE)
    #return(1)
}
    
get_new_generation <- function(i,current_pop,current_spec_ind,current_geography,migration_matrix,colors,ploidy,new_parameters,plots=F){
    lengthCurrentPop <- length(current_pop)
    if (ploidy==1) {
        #Individuals will have at least 1 offpsring plus a random poisson number.
        #n_descendents = rep(1,lengthCurrentPop) + rpois(lengthCurrentPop,growth_rate-1) 
        #n_descendents = rep(1,lengthCurrentPop) + floor(rpois(lengthCurrentPop,growth_rate-1)*(1-lengthCurrentPop/1000)) #with carrying capacity
        #growth_rate_vec <- assignParameters2species(new_parameters[[3]],current_spec_ind)
        growth_rate_vec <- new_parameters[[3]][current_spec_ind]
        #n_descendents = rep(1,lengthCurrentPop) + floor(rpois(lengthCurrentPop,growth_rate_vec-1)) #
        n_descendents = floor(rpois(lengthCurrentPop,growth_rate_vec))
        n_offspring = sum(n_descendents) # grown pop
        #The phenotypes of the new generation will be taken from a normal distribution.
        #new_pop= abs( rnorm(n_offspring,mean=rep(current_pop,n_descendents), sd = sqrt(sig2)) )
        #sig2_vec <- assignParameters2species(new_parameters[[4]],current_spec_ind)
        sig2_vec <- new_parameters[[4]][current_spec_ind]
        new_pop <- rnorm(n_offspring,mean=rep(current_pop,n_descendents), sd = rep(sqrt(sig2_vec),n_descendents))
        #new_pop= rnorm(n_offspring,mean=rep(current_pop,n_descendents), sd = sqrt(sig2))
        # reflect back values exceeding max limit
        #new_pop[new_pop>max_trait_val] = max_trait_val-abs(max_trait_val-new_pop[new_pop>max_trait_val])	
        new_pop_spec_ind = rep(current_spec_ind,n_descendents)
        new_geography <- rep(current_geography,n_descendents) #geography after reproduction, before migration.
    } else if (ploidy==2) {
        new_pop <- c()
        new_pop_spec_ind <- c()
        new_geography <- c()
        #For the diploid case we loop over each single species j.
        for (j in unique(current_spec_ind)) {
            indexes_sp_j <- which(current_spec_ind==j) #We select the indices corresponding to species j.
            current_pop_j <- current_pop[indexes_sp_j]
            current_geography_j <- current_geography[indexes_sp_j]
            couples <- gen_couples(current_pop_j,current_geography_j)
            moms <- couples[[1]]
            dads <- couples[[2]]
            minl <- couples[[3]] #the minimum length btw moms and dads determine the number of possible couples.
            couplesPositions_j <- couples[[4]]
            if (minl==0) {
                #goes extinct
                next
            } else {
                offspring_j <- get_joint_mu_s2(moms,sig2,dads,sig2)
                lengthOffspring_j <- minl
                #The following steps follow a similar pattern than in the haploid case,
                #once the number of couples "minl" is stablished.
                #Each couple will have at least 2 offspring plus a random poisson number.
                n_descendents_j <- rep(2,lengthOffspring_j) + rpois(lengthOffspring_j,2*(growth_rate-1)) #The rate is set so that we keep the user-given birth rate.
                n_offspring_j <- sum(n_descendents_j) # grown pop
                #The phenotypes of the new generation will be taken from a normal distribution.
                new_pop <- c(new_pop,abs( rnorm(n_offspring_j,mean=rep(offspring_j[[1]],n_descendents_j), sd = sqrt(offspring_j[[2]][1])) )) #why abs?
                #new_pop <- c(new_pop,rnorm(n_offspring_j,mean=rep(offspring_j[[1]],n_descendents_j), sd = sqrt(offspring_j[[2]][1])) )
                new_pop_spec_ind <- c(new_pop_spec_ind,rep(j,n_offspring_j))
                #new_geography <- c(new_geography,seq(from=1+max(current_geography_j),length.out=n_offspring_j))
                new_geography_j <- rep(couplesPositions_j,n_descendents_j)
                #migrants <- c(rep(0,length(new_geography)),sample(c(0,1,-1),size=length(new_geography_j),replace=TRUE,prob=c((1-2*probMigration),probMigration,probMigration)))
                new_geography <- c(new_geography,new_geography_j) #geography after reproduction, before migration.
                #new_geography <- new_geography+migrants
                #print(c(n_offspring_j,2*minl,lMoms+lDads,n_offspring_j/(2*minl)))
            }
        }
        
    }
    
    #Geography, migration
    P_migration_species_specific <- new_parameters[[5]]

    P_migration_species_specific_individuals <- P_migration_species_specific[new_pop_spec_ind]
    Prob_moving_individual <- Prob_moving[new_geography]
    Final_prob_migration <- Prob_moving_individual*P_migration_species_specific_individuals
    indices_inds_migrating <- which(runif(length(Final_prob_migration))<Final_prob_migration)
    areas <- 1:dim(migration_matrix)[1]
    backup_new_geography <- new_geography
    if (length(indices_inds_migrating)>0) {
        #print(paste(length(indices_inds_migrating),"individual(s) migrating"))
        for (x in 1:length(indices_inds_migrating)) new_geography[indices_inds_migrating[x]]=sample(areas,1,p=migration_matrix_transformed[new_geography[indices_inds_migrating[x]],])
    }
    #new_geography <- geography(new_geography,new_pop_spec_ind,probMigration,migration_matrix,P_migration_species_specific,typeGeography)

    # calculate death probability based on distance from other individuals
    lenNewPop <- length(new_pop)
    #death_prob_by_distance = c()
    

    #Big lengthy loop was here before
        
    #new_pop_mat <- rbind(new_pop,new_pop_spec_ind)
    #print("Testing apply ..")
    #death_prob_by_distance <- apply(new_pop_mat,2,calc_death_rate_2,new_pop,new_pop_spec_ind)
    #print("Tested!")

    if (landscape==0) {
        death_prob_by_fitness = rep(min_death_probability_by_fitness/100,lenNewPop)
    #} else if (landscape==1) {
    #    d <- dnorm(new_pop,mean=trait_optimum,sd=sd_optimum)
    #    a <- dnorm(trait_optimum,mean=trait_optimum,sd=sd_optimum)
    #    death_prob_by_fitness =  1 - (d/a)
    } else if (landscape==1) {
        d <- dnorm(new_pop,mean=trait_optimum[new_geography],sd=sd_optimum[new_geography])
        a <- dnorm(trait_optimum[new_geography],mean=trait_optimum[new_geography],sd=sd_optimum[new_geography])
        death_prob_by_fitness =  1 - (d/a)
    }

    #list_traitValues_unsorted=  get_species_trait_values_unsorted(new_pop,new_pop_spec_ind)
    #sp_means = list_traitValues_unsorted[[1]]
    #sp_index = list_traitValues_unsorted[[3]]
    #temp     = table(new_pop_spec_ind) #list_traitValues_unsorted[[4]]
    temp <- table(factor(new_pop_spec_ind, levels = 1:new_species_id))
    index_existing_species <- as.numeric(names(temp))
    #print(paste("index_exisiting_species",paste(index_existing_species,collapse=",")))
    #min_death_probability_by_fitness_intraSP = min_death_probability_by_fitness
    baseline_death_probability <- min_death_probability_by_fitness

    mean_norms <- numeric(length(unique(new_pop_spec_ind)))
    for (sp_index in unique(new_pop_spec_ind)) {
        mean_norms[sp_index] <- mean(new_pop[new_pop_spec_ind==sp_index])
    }
    
    #mean_norms = sp_means[list_traitValues_unsorted[[5]]] #
    d <- dnorm(new_pop,mean=mean_norms[new_pop_spec_ind],sd=new_parameters[[2]][new_pop_spec_ind]) #new_parameters[[2]] corresponds to the intraSP_sd_optimum per species.
    a <- dnorm(mean_norms[new_pop_spec_ind],mean=mean_norms[new_pop_spec_ind],sd=new_parameters[[2]][new_pop_spec_ind])
    #death_prob_intraspecific <- (1 - (d/a * (1-min_death_probability_by_fitness_intraSP)))
    death_prob_intraspecific <- (1 - (d/a))

    
    #print(c(list_traitValues_unsorted[[4]],table(new_pop_spec_ind)))
    
    #regressionTest <- lm(fitness_density(new_pop)~seq(1,lenNewPop))
    #slopeTest <- regressionTest$coefficients[2]
    
    #if (slopeTest<0.000001) { #If the fitness landscape is a flat line.
    #    death_prob_by_fitness = rep(min_death_probability_by_fitness,lenNewPop)
    #} else {
    #    print(paste(i,"in else!!"))
    #    for (j in 1:lenNewPop){ #THIS IS THE LOOP THAT TAKES AGES
    #        j_value = new_pop[j]
    #        j_species = new_pop_spec_ind[j]
    #        delta_j = calc_death_rate(j_value,j_species,new_pop,new_pop_spec_ind)
    #        #death_prob_by_distance = c(death_prob_by_distance,delta_j)
    #        death_prob_by_distance[j] <- delta_j
    #    }
    #   death_rate_by_fitness = log(1/exp(fitness_density(new_pop)))
    #    death_prob_by_fitness = (death_rate_by_fitness)/(max_ex_prob_fitness-min_ex_prob_fitness) * min_death_probability_by_fitness
    #}
    
    # # restrict distance to allowed vector
    # #min_distance = 0.01
    # #max_distance = 4
    # distance_vec[distance_vec<min_distance]=min_distance
    # distance_vec[distance_vec>max_distance] =max_distance
    # death_rate_by_distance =  1/distance_vec
    # death_prob_by_distance =  (death_rate_by_distance)/(1/min_distance - 1/max_distance ) * max_extinction_probability_by_distance
	
    #death_prob_Kcapacity <- (growth_rate-1)*temp/new_parameters[[1]] #death probability per species, with carryin capacity K: At carrying capacity the death rate should be equal to the growth rate.

    new_geography_factored <- factor(new_geography, levels = 1:length(K_areas))
    new_pop_spec_ind_factored <- factor(new_pop_spec_ind, levels = 1:new_species_id)
    new_geo_pop_spec <- data.frame(new_geography_factored,new_pop_spec_ind_factored)

    #free_space_area <- numeric(4)
    #for (ii in 1:4) {
    #    n_occupied <- length(subset(new_geo_pop_spec,new_geo_pop_spec$new_geography==ii)[,1])
    #    free_space_area[ii] <- K_areas[ii] - n_occupied
    #}
    temp_n_occupied <- table(new_geography_factored)
    free_space_area <- K_areas - temp_n_occupied
    free_space_area[free_space_area<0] <- 1
    Matrix_ind_per_area <- t(table(new_geo_pop_spec))
    Matrix_death_prob_K <- matrix(mapply(function(xn,species_id,area_id) (new_parameters[[3]][species_id]-1)*min(1,xn/min(new_parameters[[1]][species_id]*converted_Ks[area_id],free_space_area[area_id])), Matrix_ind_per_area, row(Matrix_ind_per_area), col(Matrix_ind_per_area)), nrow = nrow(Matrix_ind_per_area))
    
    #death_prob_Kcapacity <- myDiag(Matrix_death_prob_K[new_geo_pop_spec$new_pop_spec_ind,new_geo_pop_spec$new_geography])

    if (useCPP==1) { #This function hasn't been corrected for indices in Matrix_death_prob_K when there are extinctions.
        death_prob_Kcapacity <- assign_death_prob_K_cpp(Matrix_death_prob_K,t(data.matrix(new_geo_pop_spec)))
    } else {
        death_prob_Kcapacity <- numeric(length(new_pop_spec_ind))
        for (spp in unique(new_pop_spec_ind)) {
            index_spp <- which(new_pop_spec_ind==spp)
            #spp_row <- which(rownames(Matrix_ind_per_area)==spp)
            #print(summary(new_geography[index_spp]))
            death_prob_Kcapacity[index_spp] <- Matrix_death_prob_K[spp,new_geography[index_spp]]
        }
    }

    #print(paste("After function",i,"length",length(death_prob_Kcapacity),"growth rates",new_parameters[[3]]))
    #death_prob_Kcapacity <- vapply(new_pop_spec_ind,FUN=assign_death_prob_K,FUN.VALUE=numeric(1),new_geography,growth_rate,new_parameters[[1]],new_geo_pop_spec,free_space_area) #death probability per species, with carryin capacity K: At carrying capacity the death rate should be equal to the growth rate.
    #prob_survival <- (1-death_prob_by_fitness)*(1-death_prob_intraspecific)*(1-death_prob_Kcapacity[new_pop_spec_ind])*(1-baseline_death_probability)
    prob_survival <- (1-death_prob_by_fitness)*(1-death_prob_intraspecific)*(1-death_prob_Kcapacity)*(1-baseline_death_probability)
    #print(death_prob_Kcapacity)

    # final death probability vector 
    death_prob <- 1 - prob_survival

    r_vec = runif(length(new_pop),0,1)
    survived_new_pop_ind = which((r_vec-death_prob)>0)
    survived_new_pop = new_pop[survived_new_pop_ind]
    survived_spec_ind = new_pop_spec_ind[survived_new_pop_ind]
    survived_geography <- new_geography[survived_new_pop_ind]
    
    if (plots==T){
        writeEffectiveMigrations(i,backup_new_geography[survived_new_pop_ind],survived_geography,survived_spec_ind)
        if (landscape==1) {
            for (aG in 1:length(K_areas)) { #For n geographic areas
                #plot(x,dnorm(x,mean=trait_optimum[aG],sd=sd_optimum[aG]),type="l",main = sprintf("Area %s, N = %s, S = %s", aG,length(current_pop[survived_geography==aG]),length(unique(current_spec_ind[survived_geography==aG]))),ylim=c(0,0.5),ylab="Density")
                #points(survived_new_pop[survived_geography==aG],dnorm(survived_new_pop[survived_geography==aG],mean=trait_optimum[aG],sd=sd_optimum[aG]),pch=19,col=colors[survived_spec_ind[survived_geography==aG]])
                barplots_per_area <- table(survived_spec_ind[survived_geography==aG])
                if (length(barplots_per_area)==0) {
                    barplots_per_area <- 0
                }
                barplot(barplots_per_area,main=paste("Area",aG),col=colors[as.numeric(names(barplots_per_area))])
            }
            mtext(sprintf("Generation %s",i), outer = TRUE, cex = 1.5)
        } else {
            plot(x,(fitness_density(x)),type="l",main = sprintf("gen. %s, N = %s, S = %s", i,length(current_pop),length(unique(current_spec_ind))))
            points(survived_new_pop,fitness_density(survived_new_pop),pch=19,col=colors[survived_spec_ind])
        }
        plot(new_pop,log10(as.numeric(death_prob)),main="Death rates",ylab=expression(log[10](delta)),col=colors[new_pop_spec_ind])
        barplot(table(survived_geography),xlab="Region",ylab="Frequency",main="Geography")

    }
    #Add plots of K_a with bubbles representing K_a (R package qgraph), with arrows representing the migration rates. Color bubbles according to the sd of the landscape.

    return(list(survived_new_pop,survived_spec_ind,survived_geography,extinctionsPerArea))
}


probSpeciation <- function(current_pop,current_spec_ind,option,lambda) {
    #This function calculates the probability of speciation per species.
    #It has three options. Option 1: probability based on the distance to the mean. Option 2: probability based on the squared distance to the mean. Option 3: all species have the same probability of speciation.
    species <- unique(sort(current_spec_ind))
    numSpecies <- length(species)
    
    if(option==3) { #All species have the same probability of speciation.
        pSpeciation <- rep(lambda,numSpecies)
        return(list(pSpeciation,species))
    }

    meanDistances <- numeric(numSpecies)
    meanSquared <- numeric(numSpecies)
    k <- 1
    for (i in species) {
        vec <- current_pop[which(current_spec_ind==i)]
        if (option==1) {
            meanDistances[k] <- mean(abs(vec-mean(vec)))
        } else if (option==2) {
            meanSquared[k] <- mean((vec-mean(vec))^2)
        }
        k <- k+1
    }
    if (option==1) {
        pSpeciation <- meanDistances/sum(meanDistances)
    } else if (option==2) {
        pSpeciation <- meanSquared/sum(meanSquared)
    } 
    #The function returns two vectors. Vector 1, the probability of speciation per species. Vector 2, the species corresponding to vector 1.
    return(list(pSpeciation,species))
}

selectSpecies <- function(phenotype_distributions) {
    numSpecies <- length(phenotype_distributions)
    numAreas <- length(phenotype_distributions[[1]])
    meanMatrix <- matrix(NA,nrow=numSpecies,ncol=2)
    varMatrix <- matrix(NA,nrow=numSpecies,ncol=2)

    for (i in 1:numSpecies) {
        means <- unlist(lapply(phenotype_distributions[[i]],mean))
        vars <- unlist(lapply(phenotype_distributions[[i]],var))
        a <- which(means==min(means,na.rm=TRUE))
        b <- which(means==max(means,na.rm=TRUE))
        meanMatrix[i,1] <- means[a]
        meanMatrix[i,2] <- means[b]
        varMatrix[i,1] <- vars[a]
        varMatrix[i,2] <- vars[b]
    }
    return(which(varMatrix==min(varMatrix),arr.ind=TRUE)) #This returns two numbers: 1) the species id with the minimum variance, and 2) the island id where these individuals occur. 
}

logistic <- function(xvec,L=1,k,x0_percentage,y0=1) {
    transformationFactor <- 100
    xvec_rescaled <- xvec/(max(xvec)-min(xvec))
    xvec_rescaled <- transformationFactor*(xvec_rescaled+abs(min(xvec_rescaled)))
    x0 <- min(xvec_rescaled) + x0_percentage*(max(xvec_rescaled)-min(xvec_rescaled))
    return( (L/(1+exp(-k*(xvec_rescaled-x0)))) * y0 )
}

initialize_species_parameters <- function(new_parameters,lenNewSpecies) {
    new_parameters[[1]] <- c(new_parameters[[1]],get_K(lenNewSpecies)) 
    new_parameters[[2]] <- c(new_parameters[[2]],get_sigma2_s()) 
    new_parameters[[3]] <- c(new_parameters[[3]],get_growth_rate()) 
    new_parameters[[4]] <- c(new_parameters[[4]],get_sigma2_G())
    new_parameters[[5]] <- c(new_parameters[[5]],get_P_migration_species_specific())
    new_parameters[[6]] <- c(new_parameters[[6]],get_species_phen_range())
    return(new_parameters)
}

speciation_event_loop_3 <- function(current_pop,current_spec_ind,current_geography,lambda,new_species_id,specType,mindInd,i_external,initialMeanPhen,spp_generalistType,generalistType_list,numSpeciationsPerArea) {
    area_for_speciation <- 0 #Will return this in case there's no speciation (Area 0 does not exist).
    actualNumSpeciesSpeciating <- 0
    speciationType <- specType #Type 1 is based on distances
    sympatricORallopatric <- runif(1)
    if (specType==24) { #mixed speciation RSS and TSS
        speciationType <- sample(c(2,4),1,prob=c(probSpecType,(1-probSpecType))) #This will set the speciation type at random between types 2 and 4, with probability probSpecType of it being type 2.
    }
    speciationBasedPhen <- 0 #When set to 1 all new species die fast. Check this.
    j <- 1 #index for genealogies
    if (speciationBasedPhen==0) { #Speciation based on a probability.
        if (i_external==10){lambda=1}
        pSpec <- probSpeciation(current_pop,current_spec_ind,3,lambda)
        numSpecies <- length(pSpec[[1]]) #Gets the total number of species.
        #pSpecLambda <- lambda*numSpecies*pSpec[[1]] #Multiplies the prob of speciation per species times lambda, times the number of species.
        pSpecLambda <- pSpec[[1]] #Simpler version, without multiplying by the number of species.
        ran <- runif(numSpecies,0,1)
        speciesSpeciating <- pSpec[[2]][which((ran<pSpecLambda)==TRUE)] #Select species to speciate based on their probabilities pSpecLambda.
     
        n_speciations <- length(speciesSpeciating)
        genealogies <- matrix(0,nrow=n_speciations,ncol=2) #This matrix will save the parent-offspring relations. That is, it will keep track of which species generated which other species.

        for (i in speciesSpeciating){
            if (speciationPerArea==1) {
                if (i_external==10) {
                    area_for_speciation <- initialArea #Makes sure there are enough individuals for the first speciation in generation i_external (there might be already some individuals in another area but they might be too few).
                    sympatricORallopatric <- 0 #Makes sure the first speciation event is sympatric.
                } else {
                    area_for_speciation <- sample(unique(current_geography),1)
                }
                write(paste("Speciating species:",i,"in area",area_for_speciation,"Speciation type",speciationType,"generation",i_external),stdout())
                numSpeciationsPerArea[i] <- numSpeciationsPerArea[i]+1
                current_pop_i <- current_pop[which(current_spec_ind==i & current_geography==area_for_speciation)]
                current_spec_ind_i <- which(current_spec_ind==i & current_geography==area_for_speciation)
                #print(length(current_spec_ind_i))
            } else {
                write(paste("Speciating species:",i,"Speciation type",speciationType,"generation",i_external),stdout())
                current_pop_i <- current_pop[which(current_spec_ind==i)]
                current_spec_ind_i <- which(current_spec_ind==i)
            }
            
            length_spec_ind_i <- length(current_spec_ind_i)
            if (length_spec_ind_i < minNumIndividualsForSpeciation) {
                write(paste0("No speciation for species ",i,". Too few individuals. Exiting speciation loop."),stdout())
                numSpeciationsPerArea[i] <- numSpeciationsPerArea[i]-1
                return(list(current_spec_ind,colors,genealogies,new_species_id,new_parameters,actualNumSpeciesSpeciating,initialMeanPhen,spp_generalistType,generalistType_list,numSpeciationsPerArea,area_for_speciation))
            }
            
            if (sympatricORallopatric<Prob_sympatric()) { #Here we decide between sympatric or allopatric speciation                 
                if (speciationType==1) { #New version: Phenotypic value to split two species will be more likely to happen between values that are more far appart.
                    m=as.matrix(dist(sort(current_pop_i)))
                    distances_i <- diag(m[-nrow(m),-1]) #gets superdiagonal elements: the distances between the sorted trait values.
                    index_spec_event <- sample(c(1:length(distances_i)),1,replace=TRUE,prob=distances_i/sum(distances_i))
                    left_side <- current_spec_ind_i[1:index_spec_event]
                    right_side <- current_spec_ind_i[(index_spec_event+1):length(current_spec_ind_i)]
                } else if (speciationType==2) { #This is the original version. Take a random phenotypic value between the min and max of that species. It probably works the same as speciationType 1.
                    current_pop_i_sorted <- sort(current_pop_i)
                    spec_event = runif(1,current_pop_i_sorted[minInd],current_pop_i_sorted[length_spec_ind_i-minInd])
                    left_side  = current_spec_ind_i[which(current_pop_i < spec_event)]
                    right_side = current_spec_ind_i[which(current_pop_i > spec_event)]
                } else if (speciationType==3) { #Species assignment at random. The two new means are expected to be the same.
                    sample_size <- floor(length_spec_ind_i/2)
                    left_side <- sample(current_spec_ind_i,sample_size)
                    right_side <- setdiff(current_spec_ind_i,left_side)
                } else if (speciationType==4) { #Species assignment at random. The two new means are expected to be the same, but the sample sizes are not necessarily the same.
                    if (length_spec_ind_i>50) { #Here, I make sure there are at least 50 invidivuals in the new species.
                        sample_size <- sample(minInd:(length_spec_ind_i-minInd),1)
                    } else { #If the current speciating species has less than 50 individuals then I go back to speciationType 3.
                        sample_size <- floor(length_spec_ind_i/2)
                    }
                    left_side <- sample(current_spec_ind_i,sample_size)
                    right_side <- setdiff(current_spec_ind_i,left_side)
                } else if (speciationType==5) { #logistic
                    k <- runif(1,0,2)
                    x0_percentage <- rbeta(1,1,1)
                    prob_species_assignment <- logistic(current_pop_i,1,k,x0_percentage,1)
                    sampleSize <- ceiling(length_spec_ind_i*x0_percentage) #The following 3 lines make sure there are at least minInd individuals in each species after speciation.
                    if ( (sampleSize+minInd) >= length_spec_ind_i ) {
                        sampleSize <- length_spec_ind_i-minInd
                    }
                    sampleSize <- max(minInd,sampleSize)
                    right_side <- sample(current_spec_ind_i,sampleSize,prob=prob_species_assignment)
                    left_side <- setdiff(current_spec_ind_i,right_side)
                    #write(paste("Ancestral mean:",round(mean(current_pop_i),digits=3),"Daughter means:",round(mean(current_pop[left_side]),digits=3),",",round(mean(current_pop[right_side]),digits=3),"k =",round(k,3),"x0 =",round(x0_percentage,3)),stdout())
                }
                write(paste("Ancestral mean:",round(mean(current_pop_i),digits=3),"Daughter means:",round(mean(current_pop[left_side]),digits=3),",",round(mean(current_pop[right_side]),digits=3)),stdout())
                individuals_of_new_species <- right_side
                speciationTypeInArea <- "Sympatric"
            } else { #This is in the case of allopatry or anagenetic speciation
                print(paste("Checking:",length(unique(current_geography[which(current_spec_ind==i)])))) #checks the species is present in at least two areas
                if (length(unique(current_geography[which(current_spec_ind==i)]))>1) { #makes sure the speciating species is present in at lest two areas for allopatry to happen.
                    individuals_of_new_species <- current_spec_ind_i
                    speciationTypeInArea <- "Allopatric"
                } else {
                    write(paste0("No speciation for species ",i,". ANAGENETIC. Exiting speciation loop."),stdout())
                    return(list(current_spec_ind,colors,genealogies,new_species_id,new_parameters,actualNumSpeciesSpeciating,initialMeanPhen,spp_generalistType,generalistType_list,numSpeciationsPerArea,area_for_speciation))
                }
            }
            new_species_id <- new_species_id+1
            lenNewSpecies <- length(individuals_of_new_species)
            write(paste("Ind new sp",lenNewSpecies,";",speciationTypeInArea,"speciation"),stdout())
            #parental_species <- current_spec_ind[left_side]
            #genealogies[j,] <- c(unique(parental_species),new_species_id)
            parental_species <- i
            genealogies[j,] <- c(parental_species,new_species_id)
        
            current_spec_ind[individuals_of_new_species] <- new_species_id
            colors <- c(colors, sample(allColors,1))

            initialMeanPhen <- c(initialMeanPhen,mean(current_pop[current_spec_ind==new_species_id])) #Global variable modified here.
            generalistType <- sample(1:length(PhenotypeSpan),1)
            spp_generalistType <- c(spp_generalistType,generalistType)

            print(paste0("Parental=",i,"; spp_generalistType[i]=",spp_generalistType[i]))
            if (spp_generalistType[i]==1) {
                generalistType_list[[1]] <- c(generalistType_list[[1]],i)
            } else if (spp_generalistType[i]==2) {
                generalistType_list[[2]] <- c(generalistType_list[[2]],i)
            } else if (spp_generalistType[i]==3) {
                generalistType_list[[3]] <- c(generalistType_list[[3]],i)
            }
            

            new_parameters <- initialize_species_parameters(new_parameters,lenNewSpecies)

            write(paste("K",new_species_id,round(new_parameters[[1]][new_species_id],3),"Sigma2_s",round(new_parameters[[2]][new_species_id],3),"g",round(new_parameters[[3]][new_species_id],3),"Sigma2_G",round(new_parameters[[4]][new_species_id],3),"P_migration",round(new_parameters[[5]][new_species_id],3)),stdout())
            actualNumSpeciesSpeciating <- actualNumSpeciesSpeciating + 1 

            
            j <- j+1
        }
    } else if (speciationBasedPhen==1) { #Speciation based on differential phenotypes per region per species
        print("Deprecated option speciationBasedPhen. Check Bitbucket if you want to retrieve it.")
    }
    
    return(list(current_spec_ind,colors,genealogies,new_species_id,new_parameters,actualNumSpeciesSpeciating,initialMeanPhen,spp_generalistType,generalistType_list,numSpeciationsPerArea,area_for_speciation)) #If this is changed then also change above, there are more return functions.
}


get_species_trait_values <- function(current_pop,current_spec_ind){ #This function gets the mean and variance phenotype for each species
	species_trait_val = c()
	species_trait_var = c()
	species_index =     c()
	
	for (sp in sort(unique(current_spec_ind))){
		m = current_pop[current_spec_ind==sp]
		species_trait_val =c(species_trait_val ,mean(m))
            species_trait_var =  c(species_trait_var ,var(current_pop[current_spec_ind==sp]))
		species_index = c( species_index,sp)
	}
	return( list(species_trait_val,species_trait_var,species_index) )
}


getMoments <- function(species_trait_evolution,species_trait_evolution_var) {
    #This function is used to generate two tables: one for the mean and one for the variance of the phenotypes per time.
    #It is used later on for plotting these results.
    l <- length(species_trait_evolution)
    dimX <- l/2 #This should be time.
    dimY <- max(unlist(species_trait_evolution[seq(2,l,2)])) #This is the number of species.
    dataMean <- matrix(NA,ncol=dimX,nrow=dimY)
    dataVar <- matrix(NA,ncol=dimX,nrow=dimY)
    traitIndices <- seq(1,l,2)
    for (i in 1:dimX) {
        k <- 1
        for ( j in species_trait_evolution[[traitIndices[i]+1]] ) {
            dataMean[j,i] <- species_trait_evolution[[traitIndices[i]]][k]
            dataVar[j,i] <- species_trait_evolution_var[[traitIndices[i]]][k]
            k <- k+1
        }
    }
    return( list(dataMean,dataVar) )
}

#Functions to generate the tree
groupParentheses <- function(vec) {
    gcommas <- paste(vec,collapse=",")
    clade <- paste0("(",gcommas,")")
    return(clade)
}

groupParenthesesLen <- function(vec,len) {
    gcommas <- paste(vec,collapse=",")
    clade <- paste0("(",gcommas,"):",len)
    return(clade)
}

genNewickBLengths <- function(allGenealogies,lifespans) {
    a <- do.call(rbind,allGenealogies)
    nonZeroRows <- apply(a, 1, function(row) all(row !=0 ))
    a <- a[nonZeroRows,]
    l <- length(a[,1])
    len <- lifespans[2:(l+1),1]

    treeNewick <- groupParenthesesLen(a[1,],len[1]-1) #at root
    #print(treeNewick)
    #assign lengths for internal nodes
    for (i in 2:l) {
        i1 <- which(a[1:i,1]==a[i,1])
        i2 <- which(a[1:i,2]==a[i,1])
        i3 <- which(a[1:i,1]==a[i,2])
        i4 <- which(a[1:i,2]==a[i,2])
        pos <- unique(sort(c(i1,i2,i3,i4),decreasing=TRUE))[2]
        replacement <- groupParenthesesLen(a[i,],(len[i]-len[pos]))
        treeNewick <- gsub(paste0("([(),])",a[i,1],"([,)])"),paste0("\\1",replacement,"\\2"),treeNewick)
        #print(treeNewick)
    }
    #assign tip lengths
    tipLengths <- numeric(l+1)

    for (i in l:1) {
        lastBranchingPos <- i+1
        for (j in l:1) {
            if (a[j,1]==a[i,2]) {
                lastBranchingPos <- a[j,2]
                break;
            }
        }
        tipLengths[a[i,2]] <- lifespans[(i+1),2]-lifespans[lastBranchingPos,1]
        #print(paste("i =",i,"lastBranchingPos =",lastBranchingPos,"length =",lifespans[(i+1),2],"-",lifespans[lastBranchingPos,1],"=",tipLengths[a[i,2]]))
    }
    #for species 1
    for (j in l:1) {
        if (a[j,1]==1) {
            lastBranchingPos <- a[j,2]
            break;
        }
    }
    tipLengths[1] <- lifespans[1,2]-lifespans[lastBranchingPos,1]
    
    for (i in 1:(l+1)) {
        treeNewick <- gsub(paste0("([(),])",i,"([,)])"),paste0("\\1",i,":",tipLengths[i],"\\2"),treeNewick)
    }
    #print(tipLengths)
    treeNewick <- paste0(treeNewick,";")
    return(treeNewick)
}

update_traitOptimum <- function(trait_optimum) {
    increment <- 0.000005 #An increment of 0.000002 per generation will move an extreme optimum of e.g. 5 up to 7 in 100,000 generations.
    sds <- trait_optimum_all[1,]/1000
    #sds <- 0
    trait_optimum <- rnorm(length(trait_optimum),mean=trait_optimum^(1+increment),sd=sds)
    return(trait_optimum)
}

write_generalistType <- function(a,seed) { #a is a list with the generalistType list
    maximum <- max(unlist(lapply(a,length))) + 1
    for (i in 1:length(a)) {
        if (i==1) {
            write(paste("Type",i),file=paste0("output_generalistType_list_",seed,".txt"))
            write(a[[i]],file=paste0("output_generalistType_list_",seed,".txt"),ncolumns=maximum,append=TRUE)
        } else {
            write(paste("Type",i),file=paste0("output_generalistType_list_",seed,".txt"),append=TRUE)
            write(a[[i]],file=paste0("output_generalistType_list_",seed,".txt"),ncolumns=maximum,append=TRUE)
        }
    }
}


######-MAIN-######

###  SIMULATION PARAMETERS  ###
option_list <- list( 
    make_option(c("-v", "--verbose"), action="store_true", default=TRUE,
                help="Print extra output [default]."),
    make_option(c("-q", "--quietly"), action="store_false", 
                dest="verbose", help="Print little output."),
    make_option(c("-o", "--output"), type="double", default=0, 
                help="Print full output [default %default]."),
    make_option(c("-r", "--rseed"), type="double", default=0, 
                help="Seed [default %default]."),
    make_option(c("-f", "--landscape"), type="double", default=1, 
                help="Fitness landscape: 0=uniform, 1=normal [default %default]."),
    make_option(c("-b", "--traitOptimum"), type="double", default=5, 
                help="Trait optimum [default %default]."),
    make_option(c("-c", "--sdOptimum"), type="double", default=runif(1,10,100), #Originally it was set to 25. 
                help="SD optimum [default %default]."),
    make_option(c("-d", "--initialValue"), type="double", default=0, 
                help="Initial value of the phenotype [default %default]."),
    make_option(c("-a", "--specType"), type="double", default=4, 
                help="Speciation type: 2:vertical random (TSS), 3:horizontal symmetric, 4:horizontal random (RSS), 5:logistic, 24:mixed [default %default]."),
    make_option(c("-y", "--probSpecType"), type="double", default=0.5, 
                help="Probability of TSS compared to RSS [default %default]."),
    make_option(c("-l", "--lambda"), type="double", default=0.0002, 
                help="Speciation rate [default %default]."),
    make_option(c("-e", "--eta"), type="double", default=0, #Suggestion: for diploid mode change this value to 0.15, for haploid to 0.15.
                help="Intra-specific competition [default %default]."),
    make_option(c("-p", "--epsilon"), type="double", default=0, #Suggestion: for diploid mode change this value to 0.0001, for haploid to 0.005.
                help="Inter-specific competition [default %default]."),
    make_option(c("-s", "--sig2"), type="double", default=0.01, 
                help="Variance parent-offspring [default %default]."),
    make_option(c("-g", "--growth_rate"), type="double", default=1.2, 
                help="Birth rate [default %default]."),
    make_option(c("-n", "--ploidy"), type="double", default=1, 
                help="Ploidy: 1:haploid, 2:diploid [default %default]."),
    make_option(c("-i", "--startingPopSize"), type="double", default=500, #Make sure this value is way lower than the minimum K.
                help="Starting population size [default %default]."),
    make_option(c("-t", "--numGenerations"), type="double", default=100000, 
                help="Number of generations to run [default %default]."),
    make_option(c("-k", "--Kcapacity"), type="double", default=2000, 
                help="Initial carrying capacity [default %default]."),
    make_option(c("-z", "--SDoptimumIntra"), type="double", default=1, 
                help="Intra-species SD optimum [default %default]."),
    make_option(c("-m", "--minNumIndividuals"), type="double", default=10, 
                help="Minimum number of individuals before speciation [default %default]."),
    make_option(c("-u", "--minInd"), type="double", default=10, 
                help="Minimum number of individuals for speciation [default %default]."),
    make_option(c("-x", "--probMigration"), type="double", default=0.0001, 
                help="Probability of migration to a different region [default %default]."),
    make_option(c("-j", "--typeGeography"), type="integer", default=1, 
                help="Geography type: 1=stepping stone; 2=island model [default %default]."),
    make_option(c("--optima_file"), type="character", default=NULL,
                help="trait optima file"),
    make_option(c("--migration_file"), type="character", default=NULL,
                help="Migration matrix file"),
    make_option(c("--K_file"), type="character", default=NULL,
                help="File with K per generation per area"),
    make_option(c("-w", "--numIslands"), type="integer", default=10, 
                help="Number of islands for the Island model scenario [default %default].")
)
opt <- parse_args(OptionParser(option_list=option_list))
output <- opt$output #full output if -o 1
specType <- opt$specType
probSpecType <- opt$probSpecType
rseed <- opt$rseed
landscape <- opt$landscape
if (!is.null(opt$optima_file)) {
    trait_optima_table <- read.csv(opt$optima_file,h=F)
    trait_optima_table <- apply(trait_optima_table,2,as.numeric)
}
if (!is.null(opt$K_file)) {
    K_table <- read.csv(opt$K_file,h=F)
    K_table <- apply(K_table,2,as.numeric)
}
if (landscape==1) {
    landscape_table <- read.csv("landscapes_per_area_2.csv",header=TRUE)
    landscape_types <- landscape_table$landscapeType
    trait_optimum <- landscape_table$mean
    sd_optimum <- landscape_table$sd
    K_areas <- landscape_table$Ka
    converted_Ks <- convert_Ks(K_areas)
    if (!is.null(opt$K_file)) {
        K_areas <- K_table[1,]
        converted_Ks <- convert_Ks(K_table[1,])
    }
} else {
    trait_optimum <- opt$traitOptimum
    sd_optimum <- opt$sdOptimum
}
initialValue <- opt$initialValue
growth_rate <- opt$growth_rate # birth rate
sig2 <- opt$sig2  # variance parent-offspring
epsilon <- opt$epsilon # inter-specific competition
eta <- opt$eta # intra-specific competition
lambda <- opt$lambda # "speciation rate"
ploidy <- opt$ploidy # ploidy
startingPopSize <- opt$startingPopSize
numGenerations <- opt$numGenerations
minNumIndividualsForSpeciation <- opt$minNumIndividuals
minInd <- opt$minInd
Kcapacity  <- opt$Kcapacity
intraSP_sd_optimum <- opt$SDoptimumIntra
probMigration <- opt$probMigration
typeGeography <- opt$typeGeography
if (!is.null(opt$migration_file)) {
    migration_table <- read.csv(opt$migration_file,h=F)
    typeGeography <- 3
}
numIslands <- opt$numIslands #for typeGeography=2 (Island model)

if (rseed==0) {
    eff_seed <- sample(1:2^15, 1)
} else {
    eff_seed <- rseed
}

##Internal test.
#landscape <- 1
#trait_optima_table <- read.csv("/Users/pduchenb/Documents/Lausanne/Project_Simulator/traitOptima_baseline.csv",header=FALSE)
#trait_optima_table <- apply(trait_optima_table,2,as.numeric)
#migration_table <- read.csv("/Users/pduchenb/Documents/Lausanne/Project_Simulator/migMat_ex.csv",h=F)
#K_table <- read.csv("/Users/pduchenb/Documents/Lausanne/Project_Simulator/K_file_2.csv",header=FALSE)
#K_table <- apply(K_table,2,as.numeric)
#K_areas <- K_table[1,]
#converted_Ks <- convert_Ks(K_table[1,])
#numGenerations <- length(trait_optima_table[,1])
#numGenerations <-10000
#eff_seed <- 16827
#lambda <- 0.0005
#typeGeography <- 3

set.seed(eff_seed)

#Initial varying parameters per species.
new_parameters <- list(NULL,NULL,NULL,NULL,NULL,NULL)
new_parameters <- initialize_species_parameters(new_parameters,startingPopSize)

PhenotypeSpan <- c(0.5,2.5,6) #These numbers represent half of the phenotypic trait span: small numbers are for specialists and larger numbers for generalists.
spp_generalistType <- 3 #Vector of species with values 1 through 3 describing the generalist type. Species 1 will always be a generalist (the value 3 corresponds to the third entry in the vector PhenotypeSpan)
generalistType_list <- rep(list(c()),3) #This list has three entries, one per generalist type. Each entry is a vector of parental species that generated the generalist type of the entry.
generalistType_list[[3]] <- 1 #The first parental species is of generalist type 3.
extinctionsPerArea <- numeric(length(K_areas))
extinctionsPerArea_duePhenRanges <- numeric(length(K_areas))

write("\n",stdout())
write(sprintf("Seed for session: %s", eff_seed),stdout())
#set.seed(eff_seed)

#starting_pop <- rnorm(startingPopSize,initialValue,0.01) # initial population
starting_pop <- rnorm(startingPopSize,trait_optimum[initialArea],0.01) # initial population


if ( opt$verbose ) {
    write("\n", stdout())
    write("INPUT PARAMETERS:", stdout())
    write(paste("Landscape =",landscape), stdout())
    write(paste("Speciation type =",specType), stdout())
    write(paste("Lambda =",lambda), stdout())
    #write(paste("Sigma2_f =",round(sd_optimum,3)), stdout())
    #write(paste("Epsilon =",epsilon), stdout())
    write(paste("Type Geography =",typeGeography), stdout())
    write(paste("Initial Sigma2_G=",sig2), stdout())
    write(paste("Initial Growth rate =",growth_rate), stdout())
    write(paste("Ploidy =",ploidy), stdout())
    write(paste("Starting population size =",startingPopSize), stdout())
    write(paste("Number of generations =",numGenerations), stdout())
    write(paste("Initial carrying capacity =",Kcapacity), stdout())
    write(paste("Mininum number of individuals before speciation =",minNumIndividualsForSpeciation), stdout())
    write(paste("Mininum number of individuals for speciation =",minInd), stdout())
    write("\n", stdout())
}

pdf(file=paste0("output_plots_a_",specType,"_f_",landscape,"_",eff_seed,".pdf"),height=10,width=10) #This file will save all plots

max_trait_val = 5000
x = seq(-20,20,length.out=10000)
#plot(x,fitness_density(x),type="l")
#points(starting_pop,fitness_density(starting_pop),pch=19,col=allColors[1])
spec_ind = rep(1,length(starting_pop))
current_pop = starting_pop
current_spec_ind = spec_ind
#current_geography <- sample(seq(from=1+10^(ceiling(log10(length(starting_pop)))),length.out=length(starting_pop),by=1))
current_geography <- rep(initialArea,length(starting_pop))
species_areas_of_origin <- initialArea

len4out <- length(current_pop)

outfile_geography <- paste0("output_geography_",eff_seed,".log")
outfile_pop <- paste0("output_pop_",eff_seed,".log")
outfile_spec_ind <- paste0("output_spec_ind_",eff_seed,".log")
outfile_species <- paste0("output_species_",eff_seed,".log")
outfile_species_diversity <- paste0("output_spDiv_",eff_seed,".log")
outfile_migration <- paste0("output_migration_",eff_seed,".log")

write("\\",outfile_geography,append=FALSE)
write("Generation 0",outfile_geography,append=TRUE)
write(current_geography,outfile_geography,ncolumns=len4out,append=TRUE)

write("\\",outfile_pop,append=FALSE)
write("Generation 0",outfile_pop,append=TRUE)
write(current_pop,outfile_pop,ncolumns=len4out,append=TRUE)

write("\\",outfile_spec_ind,append=FALSE)
write("Generation 0",outfile_spec_ind,append=TRUE)
write(current_spec_ind,outfile_spec_ind,ncolumns=len4out,append=TRUE)

write("\\",outfile_migration,append=FALSE)
write("Generation 0",outfile_migration,append=TRUE)
write("0",outfile_migration,append=TRUE)

if (output==1) { #This will start writing the full information to the output file.
    outfile <- paste0("output_",eff_seed,".log")
    write("\\",outfile,append=FALSE)
    write("Generation 0",outfile,append=TRUE)
    write(current_pop,outfile,ncolumns=len4out,append=TRUE)
    write(current_spec_ind,outfile,ncolumns=len4out,append=TRUE)
}

#Tree
allGenealogies <- list()

# lowest possible fitness
#max_ex_prob_fitness = max(1/exp(fitness_density(x)))
#min_ex_prob_fitness = min(1/exp(fitness_density(x)))
min_death_probability_by_fitness = 0.00001 # lowest prob of extinction determined by fitness if it's a curve. If the fitness landscape is flat then this is the minimum death probability (when you are at maximum fitness).

pop_size=c(length(current_pop))
pop_size_species <- matrix(0,nrow=1,ncol=numGenerations) #initialize species population size matrix 
sp_div <- c(1)

list_traitValues = get_species_trait_values(current_pop,current_spec_ind) #Here we obtain the first set of species mean and variance.
species_trait_evolution = c(list(list_traitValues[[1]],list_traitValues[[3]]))
species_trait_evolution_var = c(list(list_traitValues[[2]],list_traitValues[[3]]))

k <- 1 #counter for the speciation events
new_species_id <- 1
generationSpeciation <- c() #initialize vector to save branch lengths
numSpeciationsPerArea <- numeric(length(K_areas))
numExtinctions_perArea <- numeric(length(K_areas))
totSumPerArea <- numeric(length(K_areas))


#Migration matrix
#islands <- seq(-ceiling(numIslands/2),ceiling(numIslands/2))
#migration_matrix <- matrix(0,nrow=length(islands),ncol=length(islands))
#colnames(migration_matrix) <- islands
#rownames(migration_matrix) <- islands

migration_matrix <- migrationMat(numIslands,probMigration,typeGeography,migration_table)
Prob_moving <- 1-diag(migration_matrix) #Global variable, calculated only once.
migration_matrix_transformed <- migration_matrix #do it once. This transformed matrix is only applied to individuals that will migrate.
diag(migration_matrix_transformed) <- 0 #do it once
migration_matrix_transformed <- migration_matrix_transformed/apply(migration_matrix_transformed,1,sum) #do it once
print_generations <- 100

if (landscape==1) {
    par(oma=c(0,0,2,0))
    layout(matrix(c(1:10,11,11,11,12,12),nrow=3,byrow=T)) #Allows for one row of plots per page, with the first panel set for fitness plots of 10 geographic areas.
    #par(mfrow=c(2,4),oma=c(0,0,2,0)) #Allows for two rows of plots per page (1 page per generation), with the first row set for fitness plots of 4 geographic areas.
} else {
    par(mfrow=c(2,4)) #for fitness plots
}
controlLastGeneration <- 0
#migration_anticounter <- c(0)
numberSpeciationEvents <- 1
initialMeanPhen <- c(mean(current_pop)) #initial vector of species initial mean phenotypes. This global variable is modified in the speciation function.
trait_optimum_all <- c()
extinctions_table <- c() #Table to save the species id and the area of extinction of that species (if it goes extinct).
#matrix_species_areas <- matrix(0,nrow=length(K_areas),ncol=maxNumSpeciesAllowed)
matrix_species_areas_two <- matrix(0,nrow=numGenerations,ncol=maxNumSpeciesAllowed*length(K_areas))

#write(as.vector(sapply(1:maxNumSpeciesAllowed,rep,length(K_areas))),file=paste0("test_numInd_perArea_perSpecies_",eff_seed,".csv"),ncolumns=maxNumSpeciesAllowed*length(K_areas),append=FALSE,sep=",")
#write(rep(1:length(K_areas),maxNumSpeciesAllowed),file=paste0("test_numInd_perArea_perSpecies_",eff_seed,".csv"),ncolumns=maxNumSpeciesAllowed*length(K_areas),append=TRUE,sep=",")

write(as.vector(sapply(1:maxNumSpeciesAllowed,rep,length(K_areas))),file=paste0("test_numInd_perArea_perSpecies_2_",eff_seed,".csv"),ncolumns=maxNumSpeciesAllowed*length(K_areas),append=FALSE,sep=",")
write(rep(1:length(K_areas),maxNumSpeciesAllowed),file=paste0("test_numInd_perArea_perSpecies_2_",eff_seed,".csv"),ncolumns=maxNumSpeciesAllowed*length(K_areas),append=TRUE,sep=",")

# MAIN LOOP
for (i in 1:numGenerations) {
    #K_areas <- K_table[i,] #delete
    #converted_Ks <- convert_Ks(K_table[i,]) #delete
    if (!is.null(opt$K_file)) {
        K_areas <- K_table[i,]
        converted_Ks <- convert_Ks(K_table[i,])
    }
    last_generation <- i
    if (controlLastGeneration==1) {
        break
    }
    numberOfSpecies <- length(unique(current_spec_ind))
    trait_optimum_all <- rbind(trait_optimum_all,trait_optimum)

    if (numberOfSpecies>maxNumExtantSpeciesAllowed) {
        write(paste("Reached maximum number of species:",maxNumExtantSpeciesAllowed),stdout())
        controlLastGeneration <- 1
        lambda <- 0
    }
    
    if (i%%print_generations==0) { #Write to output every print_generations generations.
        if (opt$verbose) {
            write(paste0("Generation=",i,"; Current species:",numberOfSpecies,"; Total species: ",numberSpeciationEvents,"; Total pop size:",sum(pop_size_species[,i-1]),"; Mean:",round(mean(current_pop),4)),stdout())
        }
    }
    if (i %% plot_every_x_gen == 0 | i==1){ #Add to the plots every x generations
    #if (i > 9990 | i==1){ #Add to the plots every x generations
        do_plots=T
    }else{do_plots=F}

    # REPRODUCTION OF INDIVIDUALS
    current_spec_ind_last <- current_spec_ind #Saving current spec ind and geography as "last" to check for extinctions later on.
    current_geography_last <- current_geography
    res = get_new_generation(i,current_pop,current_spec_ind,current_geography,migration_matrix,colors,ploidy,new_parameters,do_plots)

    
    #Fill in vector with extinctions per area.
    extinctSpecies <- setdiff(unique(current_spec_ind),unique(res[[2]]))
    if (length(extinctSpecies)>0) {
        for (jj in extinctSpecies) {
            areas_with_extinctions <- unique(current_geography[which(current_spec_ind==jj)])
            extinctionsPerArea[areas_with_extinctions] <- extinctionsPerArea[areas_with_extinctions]+1
        }
        print(paste("EXTINCTION OF SPECIES:",paste(extinctSpecies,collapse=","),"- IN AREA(S):",paste(areas_with_extinctions,collapse=",")))
    }
    
    current_pop = res[[1]]
    current_spec_ind = res[[2]]
    current_geography <- res[[3]]

    #Eliminate individuals outside species phenotypic ranges
    if(length(unique(current_spec_ind))>2) {
        #print(paste("spp",unique(current_spec_ind)))
        indices_to_eliminate <- c()
        for (spp in unique(current_spec_ind)) {
            #inds_within_ranges_logical <- data.table::between(current_pop,initialMeanPhen[spp]-new_parameters[[6]][spp],initialMeanPhen[spp]+new_parameters[[6]][spp]) #gets a logical vector from the whole current_pop stating which individuals are inside or outside the specified range. The indicies specific to species spp are obtained in the next line.
            inds_within_ranges_logical <- data.table::between(current_pop,initialMeanPhen[spp]-PhenotypeSpan[spp_generalistType[spp]],initialMeanPhen[spp]+PhenotypeSpan[spp_generalistType[spp]])
            #if (runif(1)>0.5) { #Each species is assigned a phenotypic range at random
            #    inds_within_ranges_logical <- data.table::between(current_pop,-2,12) #Entire phenotypic range
            #} else {
            #    inds_within_ranges_logical <- data.table::between(current_pop,0,0.3) #Small phenotypic range
            #}
            indices_to_eliminate <- c( indices_to_eliminate,intersect(which(inds_within_ranges_logical==FALSE),which(current_spec_ind==spp)) )
            #indices_to_eliminate <- c( indices_to_eliminate,subset(inds_within_ranges_logical,inds_within_ranges_logical==FALSE) )
        }
        #print(paste("A",length(indices_to_eliminate)))
        if (length(indices_to_eliminate)>0) {
            #Fill in vector with extinctions per area due to phenotype ranges.
            extinctSpecies_duePhenRanges <- setdiff(unique(current_spec_ind),unique(current_spec_ind[-indices_to_eliminate]))
            if (length(extinctSpecies_duePhenRanges)>0) {
                print(paste("EXTINCTION IN AREA(S) (due to phenotype ranges):",paste(extinctSpecies_duePhenRanges,collapse=",")))
                for (jj in extinctSpecies_duePhenRanges) {
                    areas_with_extinctions_duePhenRanges <- unique(current_geography[which(current_spec_ind==jj)])
                    extinctionsPerArea_duePhenRanges[areas_with_extinctions_duePhenRanges] <- extinctionsPerArea_duePhenRanges[areas_with_extinctions_duePhenRanges]+1
                }
            }
            current_pop <- current_pop[-indices_to_eliminate]  
            current_spec_ind <- current_spec_ind[-indices_to_eliminate]
            current_geography <- current_geography[-indices_to_eliminate]
        }
    }

    #Check for extinctions and save the area of extinction
    area_of_extinction <- c()
    extinct_species <- setdiff(unique(current_spec_ind_last),unique(current_spec_ind))
    if (length(extinct_species)>0) {
        for (aa in 1:length(extinct_species)) area_of_extinction <- c(area_of_extinction,unique(current_geography_last[which(current_spec_ind_last==extinct_species[aa])]))
        extinctions_table <- rbind(extinctions_table,c(paste(extinct_species,collapse=","),paste(area_of_extinction,collapse=",")))
    }

    #Generate vectors with local extinctions per area
    for (ss in 1:new_species_id) {
        presencePerGeneration <- table(factor(current_geography[which(current_spec_ind==ss)],levels=1:length(K_areas)))
        presencePerGeneration[presencePerGeneration>0] <- 1
        totSumPerArea <- totSumPerArea + presencePerGeneration

        presencePerGeneration_last <- table(factor(current_geography_last[which(current_spec_ind_last==ss)],levels=1:length(K_areas)))
        presencePerGeneration_last[presencePerGeneration_last>0] <- 1
        presencePerGeneration_matrix <- apply(rbind(presencePerGeneration_last,presencePerGeneration),2,as.numeric)
        indices_local_extinctions <- which(apply(presencePerGeneration_matrix,2,identical,c(1,0))==TRUE) #Vector c(1,0) represents a local extinction in a given area
        numExtinctions_perArea[indices_local_extinctions]  <- numExtinctions_perArea[indices_local_extinctions]+1     
    }
    
    #Demographic control
    # demControl <- 0
    # if (i<0) {
    #     for (j in unique(current_spec_ind)) {
    #         indicesSpecies_j <- which(current_spec_ind==j)
    #         if( length(indicesSpecies_j) > 1000 & demControl==1 ) { #Bottlenecks.
    #             write(paste("Demographic control for species",j),stdout())
    #             indicesSpecies_j_toDie <- sample(indicesSpecies_j,500)
    #             current_pop <- current_pop[-indicesSpecies_j_toDie]
    #             current_spec_ind <- current_spec_ind[-indicesSpecies_j_toDie]
    #         }
    #     }
    # }
    
    len4out <- length(current_pop)
    #print(paste("num ind:",len4out))

    if (i%%plot_every_x_gen==0) { #This will write the geographic position, traits, and indices of each individual every 'plot_every_x_gen' generations
        write("\\",outfile_geography,append=TRUE)
        write(paste("Generation",i),outfile_geography,append=TRUE)
        write(current_geography,outfile_geography,ncolumns=len4out,append=TRUE)

        write("\\",outfile_pop,append=TRUE)
        write(paste("Generation",i),outfile_pop,append=TRUE)
        write(current_pop,outfile_pop,ncolumns=len4out,append=TRUE)

        write("\\",outfile_spec_ind,append=TRUE)
        write(paste("Generation",i),outfile_spec_ind,append=TRUE)
        write(current_spec_ind,outfile_spec_ind,ncolumns=len4out,append=TRUE)

    }
    
    if (output==1) { #This will write the full information to the output file.
        write("\\",outfile,append=TRUE)
        write(paste("Generation",i),outfile,append=TRUE)
        write(current_pop,outfile,ncolumns=len4out,append=TRUE)
        write(current_spec_ind,outfile,ncolumns=len4out,append=TRUE)
    }
    
    pop_size= c(pop_size, length(current_pop))
    
    sp_div  = c(sp_div  , length(unique(current_spec_ind)))

    #if (do_plots==T){
        #plot(sp_div,type="l",main="Species diversity",xlab="t",lwd=2)		
    #}

    list_traitValues = get_species_trait_values(current_pop,current_spec_ind) #This will get the new set of species mean and variances
    species_trait_evolution = c(species_trait_evolution,  c(list(list_traitValues[[1]],list_traitValues[[3]])))
    species_trait_evolution_var <- c(species_trait_evolution_var,  c(list(list_traitValues[[2]],list_traitValues[[3]])))


    #species_trait_evolution = c(species_trait_evolution,  get_species_trait_value(current_pop,current_spec_ind))
    #species_trait_evolution_var <- c(species_trait_evolution_var,  get_species_trait_variance(current_pop,current_spec_ind))

	
    # SPECIATION 
    #if (i==(numGenerations-3000)){lambda=0}
    res_sp = speciation_event_loop_3(current_pop,current_spec_ind,current_geography,lambda,new_species_id,specType,minInd,i,initialMeanPhen,spp_generalistType,generalistType_list,numSpeciationsPerArea)
    numberNewSpecies <- res_sp[[6]]
    numberSpeciationEvents <- numberSpeciationEvents+numberNewSpecies
    current_spec_ind = res_sp[[1]]
    colors = res_sp[[2]]
    allGenealogies[[k]] <- res_sp[[3]]
    generationSpeciation <- c(generationSpeciation,i)
    k <- k+1
    new_species_id <- res_sp[[4]]
    new_parameters <- res_sp[[5]]
    initialMeanPhen <- res_sp[[7]]
    spp_generalistType <- res_sp[[8]]
    generalistType_list <- res_sp[[9]]
    numSpeciationsPerArea <- res_sp[[10]]
    if (numberNewSpecies>0) {
        species_areas_of_origin <- c(species_areas_of_origin,rep(0,numberNewSpecies)) 
        species_areas_of_origin[new_species_id] <- res_sp[[11]]
    }
    pop_size_species <- rbind(pop_size_species,matrix(0,nrow=numberNewSpecies,ncol=numGenerations)) #Increase the length of the population size matrix depending on the number of new species generated
    
    if (i>1 & sum(pop_size_species[,i-1])==0) {#If all populations go extinct.
        print("ALL SPECIES ARE EXTINCT!")
        break
    }

    for (spec_ind_forSize in unique(current_spec_ind)) { #loop to fill in population size matrix for each species
        df <- as.data.frame(table(factor(current_spec_ind, levels = 1:new_species_id)))
        pop_size_species[spec_ind_forSize,i] <- df[which(df[,1]==spec_ind_forSize),2]

        #matrix_species_areas[,spec_ind_forSize] <- table(factor(current_geography[which(current_spec_ind==spec_ind_forSize)],levels=1:length(K_areas)))
        
        placesInMatrix <- ((spec_ind_forSize-1)*length(K_areas)+1):((spec_ind_forSize-1)*length(K_areas)+length(K_areas))
        matrix_species_areas_two[i,placesInMatrix] <- table(factor(current_geography[which(current_spec_ind==spec_ind_forSize)],levels=1:length(K_areas)))
    }

    #write(c(matrix_species_areas),file=paste0("test_numInd_perArea_perSpecies_",eff_seed,".csv"),ncolumns=maxNumSpeciesAllowed*length(K_areas),append=TRUE,sep=",")
    #trait_optimum <- update_traitOptimum(trait_optimum)
    trait_optimum <- trait_optima_table[i,]
    #trait_optimum <- rnorm(length(trait_optimum),mean=trait_optimum_all[1,],sd=trait_optimum_all[1,]/1000) #only to test for varying SD.
 
}
# END OF MAIN LOOP

write.table(matrix_species_areas_two,file=paste0("test_numInd_perArea_perSpecies_2_",eff_seed,".csv"),quote=FALSE,row.names=FALSE,col.names=FALSE,append=TRUE,sep=",")

write.table(numSpeciationsPerArea,file=paste0("output_numSpeciationsPerArea_",eff_seed,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(extinctionsPerArea,file=paste0("output_numExtinctionsPerArea_",eff_seed,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(extinctionsPerArea_duePhenRanges,file=paste0("output_numExtinctionsPerArea_duePhenRanges_",eff_seed,".txt"),quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(data.frame(numExtinctions_perArea,as.numeric(totSumPerArea)),file=paste0("output_tableLocalExtinctionsPerArea_",eff_seed,".txt"),quote=FALSE,row.names=FALSE,col.names=c("localExtinctions","totalPresence"))

#Print trait optima per area
write.table(trait_optimum_all,file=paste0("output_traitOptima_perArea_",eff_seed,".csv"),quote=FALSE,row.names=FALSE,col.names=FALSE)

#Finalize matrix pop_size_species
finalSpecNum <- length(pop_size_species[,1])
pop_size_species <- cbind(rep(0,finalSpecNum),pop_size_species)
pop_size_species[1,1] <- length(spec_ind)

#Get lifespans per species
numSpecies <- length(pop_size_species[,1])
firstAppereance <- numeric(numSpecies)
lastAppereance <- numeric(numSpecies)
for (j in 1:numSpecies) {
    for (i in 1:(last_generation+1)) {
        if (pop_size_species[j,i]!=0) {
            firstAppereance[j] <- i
            break
        }
    }
    for (i in (last_generation+1):1) {
        if (pop_size_species[j,i]!=0) {
            lastAppereance[j] <- i
            break
        }
    }
}
lifespans <- data.frame(firstAppereance,lastAppereance,lastAppereance-firstAppereance+1)
names(lifespans) <- c("first","last","lifespan")
write.table(lifespans,file=paste0("output_lifespans_a_",specType,"_f_",landscape,"_",eff_seed,".csv"),quote=FALSE,row.names=FALSE,col.names=FALSE)

#Get mean and var phenotypes.
data <- getMoments(species_trait_evolution,species_trait_evolution_var)
selected_generations <- seq(1,last_generation,print_generations)
speciation_points <- as.numeric(lifespans[,1])
speciation_points_plus1 <- speciation_points+1
toKeep <- unique(sort(c(selected_generations,last_generation,speciation_points,speciation_points_plus1)))

#Print reduced output
data[[1]] <- rbind(1:length(data[[1]][1,]),data[[1]])
data[[2]] <- rbind(1:length(data[[2]][1,]),data[[2]])
pop_size_species <- rbind(1:length(pop_size_species[1,]),pop_size_species) #adds a line with generation numbers.

dataRed_means <- data[[1]][,toKeep]
dataRed_vars <- data[[2]][,toKeep]
dataRed_popSizes <- pop_size_species[,toKeep]
if (output==0) {
    write.table(t(dataRed_means),file=paste0("output_means_a_",specType,"_f_",landscape,"_",eff_seed,".csv"),quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(t(dataRed_vars),file=paste0("output_vars_a_",specType,"_f_",landscape,"_",eff_seed,".csv"),quote=FALSE,row.names=FALSE,col.names=FALSE)
    write.table(t(dataRed_popSizes),file=paste0("output_popSizes_a_",specType,"_f_",landscape,"_",eff_seed,".csv"),quote=FALSE,row.names=FALSE,col.names=FALSE)
}

lenGen <- length(data[[1]][1,])
lenData <- length(data[[1]][,1])
            
if (numberSpeciationEvents>2) {    
    #Write Newick file
    treeFileAll <- file(paste0("output_newickTreeAll_a_",specType,"_f_",landscape,"_",eff_seed,".tre"))
    writeLines(genNewickBLengths(allGenealogies,lifespans), treeFileAll)
    close(treeFileAll)
    treeAll <- read.tree(paste0("output_newickTreeAll_a_",specType,"_f_",landscape,"_",eff_seed,".tre"))
    par(mfrow=c(1,1))
    plot(treeAll,root.edge=TRUE,cex=0.5)
    axisPhylo(root.time=lifespans[2,1],backward=FALSE)

    #Write tree and traits for only tip data
    traitsTips <- data[[1]][2:lenData,lenGen]
    names(traitsTips) <- seq(1:length(traitsTips))
    toDrop <- which(is.na(traitsTips)=='TRUE')
    if (length(toDrop)>0) {
        traitsTips <- traitsTips[-toDrop]
    }
    nc=name.check(treeAll,traitsTips)
    #if (nc!="OK") {
    if (length(nc)>1) {
        tree <- drop.tip(treeAll,nc$tree_not_data)
        tree <- rescale(tree,model="depth",1)
    } else {
        tree <- rescale(treeAll,model="depth",1)
    }
    write.tree(tree,file=paste0("output_newickTreeTips_a_",specType,"_f_",landscape,"_",eff_seed,".tre"))
    write.table(traitsTips,file=paste0("output_traitsTips_a_",specType,"_f_",landscape,"_",eff_seed,".txt"),row.names=TRUE,col.names=FALSE,quote=FALSE)
} else {
    print("No tree was written, number of species is <= 2")
}


#Print species specific parameters.
write_generalistType(generalistType_list,eff_seed)
titles_species_log <- c("Species","GeneralistType","K","Sigma2_s","g","Sigma2_G","P_migration","phenRangeDiff","N_0","N_mean","lifespan","areas","area_of_origin","area_of_extinction")
write(titles_species_log,file=outfile_species,ncolumns=length(titles_species_log),sep="\t")
in_area_list <- list()
species_areas_of_extinction <- as.character(numeric(length(species_areas_of_origin))) #Initializes a character vector of zeroes.
for (z in 1:length(extinctions_table[,1])) species_areas_of_extinction[as.numeric(extinctions_table[z,1])] <- extinctions_table[z,2]

for (z in 1:length(new_parameters[[1]])) {
    initial_pop_size_species_index <- which(pop_size_species[z+1,]!=0)[1]
    mean_pop_size_species <- mean(pop_size_species[(z+1),(initial_pop_size_species_index:numGenerations+1)])
    lifespan_out <- (numGenerations+1)-initial_pop_size_species_index
    #if (pop_size_species[(z+1),(last_generation+1)]==0) {
    #    extinct <- 1
    #} else {
    #    extinct <- 0
    #}
    in_area_list[[z]] <- unique(current_geography[current_spec_ind==z])
    in_area <- paste(in_area_list[[z]],collapse=",")
    write(c(z,spp_generalistType[z],round(new_parameters[[1]][z],3),round(new_parameters[[2]][z],3),round(new_parameters[[3]][z],3),
            round(new_parameters[[4]][z],3),round(new_parameters[[5]][z],10),round(new_parameters[[6]][z],3),
            pop_size_species[z+1,initial_pop_size_species_index],mean_pop_size_species,lifespan_out,in_area,species_areas_of_origin[z],species_areas_of_extinction[z]
            ),file=outfile_species,append=TRUE,ncolumns=length(titles_species_log),sep="\t")
}


par(mfrow=c(2,2))
#Plot means and variances.
for ( i in 2:lenData ) { #i starts at 2 because i=1 has the generation numbers.
    if (i==2) {
	    m = min(as.vector(unlist(data[[1]][2:lenData,])),na.rm=TRUE)
	    M = max(as.vector(unlist(data[[1]][2:lenData,])),na.rm=TRUE)
        plot( toKeep,data[[1]][i,toKeep],type="l",col=colors[i-1],xlab="t",ylab="Mean phenotype",ylim=c(m-0.05*(M-m), M+0.05*(M-m)) )
    } else {
        lines(toKeep,data[[1]][i,toKeep],col=colors[i-1]) #i-1 in colors to correct the offset at the beginning of the loop (same below).
    }
}
for ( i in 2:length(data[[2]][,1]) ) {
    if (i==2) {
        plot(toKeep,data[[2]][i,toKeep],type="l",col=colors[i-1],xlab="t",ylab="Var phenotype",ylim=c(0,1.05*max(data[[2]][2:lenData,],na.rm=TRUE)))
    } else {
        lines(toKeep,data[[2]][i,toKeep],col=colors[i-1])
    }
}

#Plot population sizes per species
for ( i in 1:length(pop_size_species[,1]) ) {
    if(i==2) {
        maxPopSize <- max(pop_size_species[2:length(pop_size_species[,1]),])
        #plot(toKeep,pop_size_species[i,toKeep],type="l",col=colors[i-1],main="Population size per species",xlab="t",ylab="Population size",ylim=c(0,max(new_parameters[[1]])+0.1*max(new_parameters[[1]])))
        plot(toKeep,pop_size_species[i,toKeep],type="l",col=colors[i-1],main="Population size per species",xlab="t",ylab="Population size",ylim=c(0,maxPopSize+0.1*maxPopSize))
        #plot(toKeep,pop_size_species[i,toKeep],type="l",col=colors[i-1],main="Population size per species",xlab="t",ylab="Population size",ylim=c(0,15000))
    } else {
        lines(toKeep,pop_size_species[i,toKeep],col=colors[i-1])
    }
}
sumPopSizes <- apply(pop_size_species[-1,],2,sum)
plot(toKeep,sumPopSizes[toKeep],type="l",main="Total population size",xlab="t",ylab="Population size")

#Plot trait distributions per species
boxplot(current_pop~current_spec_ind, outline=F, lty=1, xlab="Species", ylab="Trait value")

#Plot variance accross species
var_no_NA <- function(x){return(var(x, na.rm=T))}
plot(toKeep,apply(data[[1]][-1,toKeep],FUN=var_no_NA,2),type="l", ylab="Var across species", xlab="Generations")

#Plot species diversity
plot(sp_div,type="l",main="Species diversity",xlab="t",lwd=2)
write.table(data.frame(1:length(sp_div),sp_div),file=outfile_species_diversity,row.names=FALSE,col.names=c("Generation","numSpecies"))

#Plot LTT
if (numberSpeciationEvents>2) {
    ltt.plot(tree,log="y",main="LTT plot")
}

#Plot Specialist-Generalist stats
barplot(as.numeric(lapply(generalistType_list,length)),names=c("Specialist","Intermediate","Generalist"),main="Speciation events")
barplot(table(factor(spp_generalistType[which(pop_size_species[,numGenerations]==0)-1],levels=1:length(generalistType_list))),names=c("Specialist","Intermediate","Generalist"),main="Extinction events")

#Plot species numbers per area
num_species_per_area <- numeric(length(K_areas))
for (y in 1:length(K_areas)) {
    num_species_per_area[y] <- length(unique(current_spec_ind[which(current_geography==y)]))
}
barplot(num_species_per_area,names=1:length(K_areas),main="Species richness per area",xlab="Area",ylab="Number of species")

#Plot moving optima
allTraitOptima <- unlist(trait_optimum_all)
minTraitOptimum <- min(allTraitOptima)-1
maxTraitOptimum <- max(allTraitOptima)+1
plot(trait_optimum_all[,1],type="l",ylim=c(minTraitOptimum,maxTraitOptimum),main="Moving optima",xlab="Generation",ylab="Optimum per area")
for (y in 2:length(trait_optimum)) lines(trait_optimum_all[,y])

#Plot areas per species
par(mfrow=c(2,1),mar=c(4,4,1,1))
ppp <- which(lapply(in_area_list,length)>0)
ppp <- ppp[match(as.numeric(tree$tip.label),ppp)]
tip=as.numeric(ppp[1])
firstHalf <- floor(length(ppp)/2)
plot(rep(1,length(in_area_list[[tip]])),in_area_list[[tip]],xlim=c(1,firstHalf),ylim=c(1,length(K_areas)),pch=15,cex=1,xaxt="n",yaxt="n",xlab="",ylab="Area",col=colors[tip])
for (i in 2:firstHalf) {
    tip=as.numeric(ppp[i])
    points(rep(i,length(in_area_list[[tip]])),in_area_list[[tip]],pch=15,cex=1,col=colors[tip])
}
axis(1,1:firstHalf,ppp[1:firstHalf],cex.axis=0.4)
axis(2,1:length(K_areas),c(1:length(K_areas)),cex.axis=0.75)
abline(v=c(1:firstHalf),lwd=0.1)

tip=as.numeric(ppp[firstHalf+1])
secondHalf <- length(seq((firstHalf+1),length(ppp)))
plot(rep(1,length(in_area_list[[tip]])),in_area_list[[tip]],xlim=c(1,secondHalf),ylim=c(1,length(K_areas)),pch=15,cex=1,xaxt="n",yaxt="n",xlab="Tip",ylab="Area",col=colors[tip])
for (i in 2:secondHalf) {
    tip=as.numeric(ppp[i+firstHalf])
    points(rep(i,length(in_area_list[[tip]])),in_area_list[[tip]],pch=15,cex=1,col=colors[tip])
}
axis(1,1:secondHalf,ppp[(firstHalf+1):length(ppp)],cex.axis=0.4)
axis(2,1:length(K_areas),c(1:length(K_areas)),cex.axis=0.75)
abline(v=c(1:secondHalf),lwd=0.1)

#Plot extant timetree with areas per species
if (numberSpeciationEvents>2) {
    for (i in 1:length(tree$tip.label)) {
        tree$tip.label[i] <- sprintf("%03d",as.numeric(tree$tip.label[i]))
        presenceAreasForTree <- table(factor(sort(in_area_list[ppp[i]][[1]]),levels=1:length(K_areas)))
        tree$tip.label[i] <- paste(tree$tip.label[i],paste(presenceAreasForTree,collapse=" "),sep=" - ")
    }
    par(mfrow=c(1,1))
    plot.phylo(tree,cex=0.4)
}

##NOTES

# distance and species ID (potentially allowing hybridization)

# 3D plots, axis 1 trait, axis 2 fitness, axis 3 geographic distance.
# For that, make that the individual also moves along an axis of distance, and we draw the values from a normal distribution centered at the current position.
# For diploids, sum the normal distributions of both parents and from the resulting new normal distribution draw the offspring value. Add option for selfing.

#Add Carrying capacity and fitness landscape per geographic area (allow for multimodal landscapes too). Look at Boucher 2017.
#Add biogeographic speciation (if individuals don't mix in x generations then they speciate).
#Add a way to compare selection coefficients and sigmas between micro and macro scales.
#Track genealogies as well (for multispecies coalescent).

dev.off()
