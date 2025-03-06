# phylotraitsim: Simulations of intra- and inter-specific trait evolution
This is code to simulate individual-based phenotypic data. For details please refer to
Duchen, Pablo, et al. "On the effect of asymmetrical trait inheritance on models of trait evolution." Systematic Biology 70.2 (2021): 376-388.

### Requirements 
The simulator is written in R requires the libraries `ape`, `geiger`, and `optparse`.
The code was tested on different MacOS and Linux distributions using R 4.x.


### Usage

1) To display the full options type in a Terminal window:

`Rscript simulator0.2.R --help`

2) To run a simulation of 25000 generations, speciation rate (lambda) 0.0005,
a migration matrix given in the file migMat_ex_2.csv,
assuming the seed is 10101 and the evolutionary rate is 0.01 we run the following:

Rscript logistic.R 10101 0.01 (This would generate the optima file for the next step), and

Rscript simulator0.2.R -t 25000 --optima_file optima_trajectories_movingMidpoint_10101.csv --lambda 0.0005 --migration_file migMat_ex_2.csv --rseed 10101

### Simulations runs for the article "Environmental fluctuations trigger faster evolution and turnover"
For the study by Duchen et al. (2025) "Environmental fluctuations trigger faster evolution and turnover" we ran simulations in an HPC computing cluster with the following code.
```
#!/bin/bash
#SBATCH -J simRun                 
#SBATCH -o simRun.%j.out                             
#SBATCH -n 1        
#SBATCH -c 1                   
#SBATCH --mem 2G            
#SBATCH -t 10:00:00       
#SBATCH --array=1-300%300

module load lang/R/4.2

parameters=seeds_rates_3.txt
p_a=$(cat $parameters | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $1}')
p_b=$(cat $parameters | awk -v var=$SLURM_ARRAY_TASK_ID 'NR==var {print $2}')

Rscript logistic.R ${p_a} ${p_b}
Rscript simulator0.3.R -t 25000 --optima_file optima_trajectories_movingMidpoint_${p_a}.csv --lambda 0.0005 --migration_file migMat_ex_2.csv --rseed ${p_a}
```

The script above will use the parameters from the file `seed_rates_3.txt` and parse them to `logistic.R` and `simulator0.3.R`. These parameters (p_a and p_b) contain a seed and three different $\sigma_E$'s corresponding to three evolutionary pace scenarios. A total of 100 simulations are run per scenario, thus, the file `seed_rates_3.txt` has 300 lines.

## Description of output files generated by logistic.R
Assuming the seed used was 10101 then the output files would be the following:
#### midpoint_trajectory_10101.csv
Contains the trajectory of the midpoint at each generation.
#### optima_trajectories_movingMidpoint_10101.csv
Contains the optima trajectories for each area.

## Description of output files generated by simulator0.2.R
#### output_generalistType_list_10101.txt
Contains species id's divided in three vectors: 
Type 1 (specialist with a phenotype span of 1), 
Type 2 (intermediate with a phenotype span of 5), and 
Type 3 (generalist with a phenotype span of 12).
The species id's are sometimes repeated when they were the parent of a new species.
Thus, this file also shows the frequency at which a given species divided in two and the generalist type it belongs to.

#### output_geography_10101.log
Area index of each individual recorded every 1000 generations.

#### output_lifespans_a_4_f_1_10101.csv
Table with one row per species, and four columns. 
Column 1: Species ID.
Column 2: Generation at which a given species originated,
Column 3: Generation at which that same species went extinct, and
Colunn 4: The lifespan in generations.

#### output_means_a_4_f_1_10101.csv
Table with mean trait values per species at certain generations.
The first column represents generations, and each subsequent column represents a species.
The mean trait values per species per generation are recorded for each species at certain generation times (around speciation events).

#### output_migration_10101.log
Information about migration events recorded every 1000 generations.
The first value represents the number of migration events.
The second set of values represent the areas of origin.
The third set of values represent the areas of destination.
The fourth set of values represent the species that migrated.
For example, the output:

Generation 1000

1 4 5 3

means that, at generation 1000, there was one migration event
from area 4 to area 5, involving species 3. Likewise, the output:

Generation 8000

3 4 4 6 5 3 7 4 4 3

means that, at generation 8000, there were 3 migration events
from areas 4, 4, and 6, to areas 5, 3, and 7, respectively,
involving species 4, 4, and 3, respectively.

#### output_newickTreeAll_a_4_f_1_10101.tre
Phylogenetic tree in Newick format including all extant and extinct species.

#### output_newickTreeTips_a_4_f_1_10101.tre
Phylogenetic tree in Newick format including only extant species.

#### output_numExtinctionsPerArea_10101.txt
Number of species that went fully extinct, and the area where it happened (one line per area).

#### output_numExtinctionsPerArea_duePhenRanges_10101.txt
Number of species that went fully extinct due to exceeding its phenotypic range, and the area where it happened (one line per area).
This file should be a subset of the more complete output file output_numExtinctionsPerArea_10101.txt

#### output_numSpeciationsPerArea_10101.txt
Number of speciations per area (one line per area).

#### output_plots_a_4_f_1_10101.pdf
Various plots summarizing the simulation run.

#### output_popSizes_a_4_f_1_10101.csv
Table with total population size per species (one column per species)
at various time points (around speciation events).
The first column indicates the generation time point.

#### output_pop_10101.log
Trait values of each individual recorded every 1000 generations.

#### output_spDiv_10101.log
Number of species present at every generation.

#### output_spec_ind_10101.log
Species index of every individual recorded every 1000 generations.

#### output_species_10101.log
Summary table for each species (rows) and different statistics (columns) including:

Species: Species ID

GeneralistType: Generalist type (1: specialist, 2: intermediate, 3: generalist). For phenotypic ranges of each type see output_generalistType.

K: Maximum carrying capacity. It is set by default to 4000, unless the population size of a new species is larger than that (in which case K will be assigned to the size of the new species).

Sigma2_s: Intra-specific sd^2 optimum.

g: Growth rate.

Sigma2_G: Parent-offspring variance.

P_migration: Species specific probability of migration, independent of the migration matrix (function get_P_migration_species_specific).

phenRangeDiff: Species specific phenotypic range, assigned at random between 0 and 5 (currently not used).

N_0: Initial population size.

N_mean: Mean population size across all generations.

lifespan: Life span in generations.

areas: Areas where the species is present.

area_of_origin: Area where the species originated.

area_of_extinction: Area where the species went extinct.

#### output_tableLocalExtinctionsPerArea_10101.txt
Table with the number of local extinctions per area (one line per area), plus the total presence per area 
(i.e. the total number of generations where there were any individuals living there).

#### output_traitOptima_perArea_10101.csv
Trait optima per area (columns) recorded every generation (rows).

#### output_traitsTips_a_4_f_1_10101.txt
Table with trait values at the tips the phylogeny (only for extant species).
Column one is the species index, and column two the mean trait value.

#### output_vars_a_4_f_1_10101.csv
Table with variance trait values per species at certain generations.
The first column represents generations, and each subsequent column represents a species.
The variance trait values per species per generation are recorded for each species at certain generation times (around speciation events).

#### test_numInd_perArea_perSpecies_2_10101.csv
Same as output_popSizes but further separated by area.
The first row indicates the species, and the second row the area.
(In this example, only the first 10000 generations have been printed).

## Scripts for post-processing some output files.
### Script analize_K.R
Takes the files test_numInd_perArea_perSpecies_2, output_newickTreeAll, and landscapes_per_area_2.csv as input, and outputs the files res_species and res_analyze_K.

#### res_species_10101.csv
This files contains the species richness per area.
Assuming there are 9 areas, this file contains the number of species 
after the end of the simulation run (first 9 values), 
after 30000 generations (next 9 values),
after 15000 generations (next 9 values), and 
after 7500 generations (next 9 values).

#### res_analyze_K_10101.csv
Contains the total number of individuals per area (one area per column) recorded at every generation.

### Script analize_full_output.R
Takes the files test_numInd_perArea_perSpecies_2, output_newickTreeAll, and landscapes_per_area_2.csv as input, and outputs the file res_analyze_full.

#### res_analyze_full_10101.csv
This files contains the total number of speciations per area (first 9 values, assuming a total of 9 areas),
then the number of extinctions per area (next 9 values),
and the total time spent (species presence) in each area (final 9 values).
