# phylotraitsim
## Usage
(Assuming the seed is 10101 and the evolutionary rate is 0.01)

Rscript logistic.R 10101 0.01

Rscript simulator0.2.R -t 25000 --optima_file optima_trajectories_movingMidpoint_10101.csv --lambda 0.0005 --migration_file migMat_ex_2.csv --rseed 10101

## Description of output files generated by logistic.R
Assuming the seed used was 10101 then the output files would be the following:
### midpoint_trajectory_10101.csv
Contains the trajectory of the midpoint at each generation.
### optima_trajectories_movingMidpoint_10101.csv
Contains the optima trajectories for each area.

## Description of output files generated by simulator0.2.R
### output_generalistType_list_10101.txt
Contains indices of parental species divided in three vectors: 
Type 1 (specialist with a phenotype span of 1), 
Type 2 (intermediate with a phenotype span of 5), and 
Type 3 (generalist with a phenotype span of 12).

### output_geography_10101.log
Area index of each individual recorded every 1000 generations.

### output_lifespans_a_4_f_1_10101.csv
Table with one row per species, and three columns. 
Column 1: Generation at which a given species originated,
Column 2: Generation at which that same species went extinct, and
Colunn 3: The lifespan in generations.

### output_means_a_4_f_1_10101.csv
Table with mean trait values per species at certain generations.
The first column represents generations, and each subsequent column represents a species.
The mean trait values per species per generation are recorded for each species at certain generation times (around speciation events).

### output_migration_10101.log
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

### output_newickTreeAll_a_4_f_1_10101.tre
Phylogenetic tree in Newick format including all extant and extinct species.

### output_newickTreeTips_a_4_f_1_10101.tre
Phylogenetic tree in Newick format including only extant species.

### output_numExtinctionsPerArea_10101.txt
Number of species that went fully extinct, and the area where it happened (one line per area).

### output_numExtinctionsPerArea_duePhenRanges_10101.txt
Number of species that went fully extinct due to exceeding its phenotypic range, and the area where it happened (one line per area).
This file should be a subset of the more complete output file output_numExtinctionsPerArea_10101.txt

### output_numSpeciationsPerArea_10101.txt
Number of speciations per area (one line per area).

### output_plots_a_4_f_1_10101.pdf
Various plots summarizing the simulation run.

### output_popSizes_a_4_f_1_10101.csv
Table with total population size per species (one column per species)
at various time points (around speciation events).
The first column indicates the generation time point.

### output_pop_10101.log
Trait values of each individual recorded every 1000 generations.

### output_spDiv_10101.log
Number of species present at every generation.

### output_spec_ind_10101.log
Species index of every individual recorded every 1000 generations.

### output_species_10101.log
Summary table for each species (rows) and different statistics (columns) including: generalist type, carrying capacity (K), various rates, presence in areas, lifespan, etc.

### output_tableLocalExtinctionsPerArea_10101.txt
Table with the number of local extinctions per area (one line per area), plus the total presence per area 
(i.e. the total number of generations where there were any individuals living there).

### output_traitOptima_perArea_10101.csv
Trait optima per area (columns) recorded every generation (rows).

### output_traitsTips_a_4_f_1_10101.txt
Table with trait values at the tips the phylogeny (only for extant species).
Column one is the species index, and column two the mean trait value.

### output_vars_a_4_f_1_10101.csv
Table with variance trait values per species at certain generations.
The first column represents generations, and each subsequent column represents a species.
The variance trait values per species per generation are recorded for each species at certain generation times (around speciation events).

## Description of output files generated by analyzeK.R