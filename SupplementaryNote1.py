### Supplementary Note 1 ###
# python code to run the coalescent simulations under the demographic model shown in Supplementary Figure 4A #
# see Supplementary Section 4 for more information #

import tskit
import msprime
import numpy as np
import math
import random
import sys


# note: populations WEA* and EEA in the code correspond to populations named WEC* and EEC throughout the article #

### function to run simulations ###
generation_time = 29
def run_simulation(sequence_length, demography, seed = np.random.randint(1, 2**32 - 1, 1)):
    ts = msprime.sim_ancestry(
        demography = demography,
        sequence_length = sequence_length,
        samples = sample_list,
        recombination_rate = 1e-8,
        random_seed = seed
    )
    return ts

### split and sample times ###
#
t_AMH_sample = int(200000 / generation_time)
t_AFR1_split = int(172000 / generation_time)
t_AFR2_split = int(100000 / generation_time)
t_AFR3_split = int(65000 / generation_time)
t_OOA_split = int(60000 / generation_time)
t_BEA_split = int(57500 / generation_time)
t_WEA_EEA_split = int(46000 / generation_time)
t_WEA2_split = int(40000 / generation_time)

# this parameter controls how many years before WEC2/WEA2 did WEC/WEA leave the hub #
# N is 43000, 42000 and 41000 in the simulations whose results are shown in Supplementary figure A, B and C respectively #
t_WEA_split = int(N / generation_time)
#
t_UST_sample = int(44377 / generation_time)
t_UST_split = t_WEA_EEA_split
#
t_TIA_sample = time = int(39565 / generation_time)
t_TIA_split = t_TIA_sample + 34
t_AR_sample = int(33592 / generation_time)
t_AR_split = t_TIA_sample + 10
#
t_KOS_sample = int(38052 / generation_time)
t_KOS_split = t_KOS_sample + 34
t_SUN_sample = int(34170 / generation_time)
t_SUN_split = t_KOS_sample + 10
#


### admixture times ###
t_BEA2WEA = int(25000 / generation_time)
t_ADM2 = int(10000 / generation_time)


### population sizes ###
# Ne values from Gravel et al. 2011 #
Ne_AFR = 14474
Ne_OOA = 1861
Ne_HUB = Ne_OOA
Ne_BEA_0 = 500
Ne_WEA_0 = 1032
Ne_EEA_0 = 550


# growth rates #
# We need to work out the starting (diploid) population sizes based on
# the growth rates provided for these two populations
# Rescale Gravel et al. 2011 parameters so that the final Ne is the same #
# Ne_final_EAS = 86674.78397785178 | Ne_0 = 550 #
# Ne_final_EUR = 40915.07868289492 | Ne_0 = 1032 #
# r = math.log(Ne/Ne_0)/t #
# r_EUR = 0.00232 | r_EAS = 0.0031899999999999997 #
r_WEA = 0.00232
r_EEA = 0.00319
r_BEA = r_WEA

Ne_WEA = Ne_WEA_0 / math.exp(-r_WEA * t_WEA_EEA_split)
Ne_EEA = Ne_EEA_0 / math.exp(-r_EEA * t_WEA_EEA_split)
Ne_BEA = Ne_BEA_0 / math.exp(-r_BEA * t_BEA_split)


### define populations ###
# 

demography = msprime.Demography()
demography.add_population(name = "AMH", initial_size = 7300, default_sampling_time = t_AMH_sample)
demography.add_population(name = "AFR1", initial_size = Ne_AFR)
demography.add_population(name = "AFR2", initial_size = Ne_AFR)
demography.add_population(name = "AFR3", initial_size = Ne_AFR)
demography.add_population(name = "OOA", initial_size = Ne_OOA)
demography.add_population(name = "HUB", initial_size = Ne_OOA, default_sampling_time = t_WEA_EEA_split + 1)
demography.add_population(name = "HUB_WEA", initial_size = Ne_OOA)
demography.add_population(name = "HUB_WEA2", initial_size = Ne_OOA)
demography.add_population(name = "BEA", initial_size = Ne_BEA, growth_rate = r_WEA, default_sampling_time = t_BEA2WEA + 1)
demography.add_population(name = "EEA", initial_size = Ne_EEA, growth_rate = r_EEA)
demography.add_population(name = "WEA", initial_size = Ne_WEA, growth_rate = r_WEA, default_sampling_time = t_BEA2WEA + 1)
demography.add_population(name = "WEA2", initial_size = Ne_WEA, growth_rate = r_WEA, default_sampling_time = t_BEA2WEA + 1)
demography.add_population(name = "UST", initial_size = 600, default_sampling_time = t_UST_sample)
demography.add_population(name = "KOS", initial_size = 600, default_sampling_time = t_KOS_sample)
demography.add_population(name = "TIA", initial_size = 600, default_sampling_time = t_TIA_sample)
demography.add_population(name = "SUN", initial_size = 600, default_sampling_time = t_SUN_sample)
demography.add_population(name = "AR", initial_size = 600, default_sampling_time = t_AR_sample)
demography.add_population(name = "WEA_BEA_a", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA_BEA_b", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA_BEA_c", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA_EEA_a", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA_EEA_b", initial_size = Ne_EEA, growth_rate = r_EEA)
demography.add_population(name = "WEA_EEA_c", initial_size = Ne_EEA, growth_rate = r_EEA)
demography.add_population(name = "WEA_BEA_EEA_a", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA_BEA_EEA_b", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA_BEA_EEA_c", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA_BEA_EEA_d", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA2_BEA_a", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA2_BEA_b", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA2_BEA_c", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA2_EEA_a", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA2_EEA_b", initial_size = Ne_EEA, growth_rate = r_EEA)
demography.add_population(name = "WEA2_EEA_c", initial_size = Ne_EEA, growth_rate = r_EEA)
demography.add_population(name = "WEA2_BEA_EEA_a", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA2_BEA_EEA_b", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA2_BEA_EEA_c", initial_size = Ne_WEA, growth_rate = r_WEA)
demography.add_population(name = "WEA2_BEA_EEA_d", initial_size = Ne_WEA, growth_rate = r_WEA)


### admixture events ###

demography.add_admixture(
    t_ADM2,
    derived = "WEA_BEA_EEA_a",
    ancestral = ["WEA_BEA_a", "EEA"],
    proportions = [0.90, 0.10],
)
demography.add_admixture(
    t_ADM2,
    derived = "WEA_BEA_EEA_b",
    ancestral = ["WEA_BEA_a", "EEA"],
    proportions = [0.80, 0.20],
)
demography.add_admixture(
    t_ADM2,
    derived = "WEA_BEA_EEA_c",
    ancestral = ["WEA_BEA_b", "EEA"],
    proportions = [0.90, 0.10],
)
demography.add_admixture(
    t_ADM2,
    derived = "WEA_BEA_EEA_d",
    ancestral = ["WEA_BEA_b", "EEA"],
    proportions = [0.80, 0.20],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA_BEA_a",
    ancestral = ["WEA", "BEA"],
    proportions = [0.90, 0.10],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA_BEA_b",
    ancestral = ["WEA", "BEA"],
    proportions = [0.80, 0.20],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA_BEA_c",
    ancestral = ["WEA", "BEA"],
    proportions = [0.70, 0.30],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA_EEA_a",
    ancestral = ["WEA", "EEA"],
    proportions = [0.75, 0.25],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA_EEA_b",
    ancestral = ["WEA", "EEA"],
    proportions = [0.5, 0.5],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA_EEA_c",
    ancestral = ["WEA", "EEA"],
    proportions = [0.25, 0.75],
)

demography.add_admixture(
    t_ADM2,
    derived = "WEA2_BEA_EEA_a",
    ancestral = ["WEA2_BEA_a", "EEA"],
    proportions = [0.90, 0.10],
)
demography.add_admixture(
    t_ADM2,
    derived = "WEA2_BEA_EEA_b",
    ancestral = ["WEA2_BEA_a", "EEA"],
    proportions = [0.80, 0.20],
)
demography.add_admixture(
    t_ADM2,
    derived = "WEA2_BEA_EEA_c",
    ancestral = ["WEA2_BEA_b", "EEA"],
    proportions = [0.90, 0.10],
)
demography.add_admixture(
    t_ADM2,
    derived = "WEA2_BEA_EEA_d",
    ancestral = ["WEA2_BEA_b", "EEA"],
    proportions = [0.80, 0.20],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA2_BEA_a",
    ancestral = ["WEA2", "BEA"],
    proportions = [0.90, 0.10],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA2_BEA_b",
    ancestral = ["WEA2", "BEA"],
    proportions = [0.80, 0.20],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA2_BEA_c",
    ancestral = ["WEA2", "BEA"],
    proportions = [0.70, 0.30],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA2_EEA_a",
    ancestral = ["WEA2", "EEA"],
    proportions = [0.75, 0.25],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA2_EEA_b",
    ancestral = ["WEA2", "EEA"],
    proportions = [0.5, 0.5],
)
demography.add_admixture(
    t_BEA2WEA,
    derived = "WEA2_EEA_c",
    ancestral = ["WEA2", "EEA"],
    proportions = [0.25, 0.75],
)


### population splits ###
demography.add_mass_migration(time = t_SUN_split, source ="SUN", dest ="KOS", proportion = 1)
demography.add_mass_migration(time = t_AR_split, source ="AR", dest ="TIA", proportion = 1)
demography.add_mass_migration(time = t_KOS_split, source ="KOS", dest ="WEA", proportion = 1)
demography.add_mass_migration(time = t_TIA_split, source ="TIA", dest ="EEA", proportion = 1)
demography.add_population_split(time = t_WEA_split, derived = ["HUB_WEA2", "WEA"], ancestral ="HUB_WEA")
demography.add_population_split(time = t_WEA2_split, derived = ["WEA2"], ancestral ="HUB_WEA2")
demography.add_population_split(time = t_WEA_EEA_split, derived = ["HUB_WEA", "EEA"], ancestral ="HUB")
demography.add_population_split(time = t_UST_split, derived = ["UST"], ancestral ="HUB")
demography.add_population_split(time = t_BEA_split, derived = ["HUB", "BEA"], ancestral ="OOA")
demography.add_mass_migration(time = t_OOA_split, source ="OOA", dest ="AFR3", proportion = 1)
demography.add_mass_migration(time = t_AFR3_split, source ="AFR3", dest ="AFR2", proportion = 1)
demography.add_mass_migration(time = t_AFR2_split, source ="AFR2", dest ="AFR1", proportion = 1)
demography.add_population_split(time = t_AFR1_split, derived = ["AFR1"], ancestral ="AMH")

demography.sort_events()



### prepare sample set and list of individual names for vcf file ###
sample_list = []
individual_names = []

for pop in ["AMH", "KOS", "TIA", "UST", "SUN", "AFR3", "EEA", "WEA", "WEA2",
            "WEA_EEA_a", "WEA_EEA_b", "WEA_EEA_c", "WEA_BEA_a", "WEA_BEA_b", "WEA_BEA_c",
            "WEA_BEA_EEA_a", "WEA_BEA_EEA_b", "WEA_BEA_EEA_c", "WEA_BEA_EEA_d",
            "WEA2_EEA_a", "WEA2_EEA_b", "WEA2_EEA_c", "WEA2_BEA_a", "WEA2_BEA_b", "WEA2_BEA_c",
            "WEA2_BEA_EEA_a", "WEA2_BEA_EEA_b", "WEA2_BEA_EEA_c", "WEA2_BEA_EEA_d"]:
    sample_list.append(msprime.SampleSet(10, ploidy = 2, population = pop))
    individual_names = individual_names + [ pop + "-" + pop + "_{0}".format(i) for i in range(1, 11)]


### simulation parameters ###
mu = 1.25*10**-8
seqlen = 3*10**9


### run simulations and save genomes in vcf file ###
ts = run_simulation(seqlen, demography)
mts = msprime.sim_mutations(ts, rate = mu) # add mutations #

with open(str(sys.argv[0][:-3]) + str(sys.argv[1]) + ".vcf", "w") as vcf_file:
    mts.write_vcf(vcf_file, individual_names=individual_names)

