import argparse
import random
import subprocess
from colorama import Fore, Style


### FUNCTIONS ###

# identify derived alleles on one individual #
def list_derived_ind(FID_IID, dataset, outdir = "", het_rndlow = 0):
    sample = FID_IID.split()[1]
    print(Fore.LIGHTYELLOW_EX + Style.BRIGHT + "\nIdentifying sites with derived allele for sample " +  sample + " ...", end = " ")
    # write file with FID and IID for plink --keep #
    with open(outdir + sample + ".txt", "w") as outfile:
        outfile.write(FID_IID)
    # run plink to get stats for each snp #
    subprocess.run(["plink", "--bfile", dataset, "--keep-allele-order", "--freqx",
                    "--keep", outdir + sample + ".txt", "--out", outdir + sample], stdout = subprocess.DEVNULL)
    # identify sites with derived allele #
    list_derived = []
    with open(outdir + sample + ".frqx") as sites:
        for site in sites:
            # if hom for derived allele (= 1 in column 7) add to list #
            if site.split()[6] == "1":
                list_derived.append(site.split()[1])
            # if het (= 1 in column 6) add to list with 0.5 probability #
            # if you want to include all derived alleles found at heterozygous genotypes set het_rndlow to 1 #
            elif site.split()[5] == "1" and random.randint(het_rndlow, 1) == 1:
                list_derived.append(site.split()[1])
    # write list of sites with derived allele to file #
    print(sample, "has", len(list_derived), "sites with derived allele", Style.RESET_ALL)
    if het_rndlow == 1: het_status = "allhet_"
    else: het_status = ""
    with open(outdir + "derived_alleles_" + het_status + sample + ".snps", "w") as outfile:
        outfile.write("\n".join(list_derived))
    # delete intermediate files #
    subprocess.run(["rm", outdir + sample + ".log", outdir + sample + ".frqx", outdir + sample + ".nosex", outdir + sample + ".txt"])
#


###

# initialize the parser #
parser = argparse.ArgumentParser(
    description = "Identify sites with derived allele - You must first create a plink bynary file where the ancestral allele is set to A1 (e.g. by using plink --reference-allele)"
)

# add arguments #
parser.add_argument("-d", "--dataset", required = True, metavar = "[file prefix]",
                    help = "prefix of input plink files")
group = parser.add_mutually_exclusive_group(required = True)
group.add_argument("-i", "--indlist", metavar = "[file]",
                   help = "list of individuals for which derived alleles have to be identified, "
                          "one per line with FID and IID of plink .fam file "
                          "(tip: add Ancestor.REF Ancestor.REF as a safe check)")
group.add_argument("-p", "--poplist", metavar = "[file]",
                   help = "list of populations for which mean number of derived alleles needs to be computed, "
                          "one per line and equal to FID of plink .fam file "
                          "(tip: add Ancestor.REF as a safe check)")
parser.add_argument("-a", "--all", metavar = "", action = "store_const", const = 1, default = 0, dest = "rndlow",
                    help = "keep all derived sites in diploid individuals"
                           "[default behaviour: keep derived sites at heterozygous sites with 0.5 probability]")

# parse the arguments #
args = parser.parse_args()

###

# if poplist is provided indentify derived alleles for every individual of each of the populations provided #
# create a directory for each population to store individual results #
if args.poplist != None:
    with open(args.poplist) as poplist:
        for pop in poplist:
            pop = pop.splitlines()[0]
            print(Fore.RED + Style.BRIGHT + "\n\nAnalysing ", pop, ":", Style.RESET_ALL)
            subprocess.run(["mkdir", pop])
            # get list of individuals for each pop #
            out = subprocess.run(["grep", "-w", "^" + pop, args.dataset + ".fam"], capture_output = True, text = True)
            # write indlist for list_derived_ind() function #
            with open(pop + "/indlist.txt", "w") as outfile:
                outfile.write(out.stdout)
            # get list of derived allele for each ind in the population #
            with open(pop + "/indlist.txt") as indlist:
                for ind in indlist:
                    list_derived_ind(ind, args.dataset, pop + "/", args.rndlow)

# if indlist is provided identify derived allele for each individual #
else:
    with open(args.indlist) as indlist:
        for ind in indlist:
            list_derived_ind(ind, args.dataset, het_rndlow = args.rndlow)
#
