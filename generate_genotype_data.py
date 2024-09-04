import itertools
from collections import defaultdict
import random

def generate_genotype_data(loci,
                           alleles,
                           Nm = None,
                           Nf = None,
                           Nm_avg = None,
                           Nm_std = None,
                           Nf_avg = None,
                           Nf_std = None,
                           covariance_avg=None,
                           covariance_std=None,
                           seed=None):
    """
    Generate all possible genotypes given loci and alleles, with random fitness values,
    and specified averages and standard deviations for the number of males and females.

    Parameters:
    loci (int): The number of loci.
    alleles (int): The number of alleles per locus.
    Nm (int, optional): The total number of males in the population. If specified, Nm_avg and Nm_std are ignored.
    Nf (int, optional): The total number of females in the population. If specified, Nf_avg and Nf_std are ignored.
    Nm_avg (int, optional): The average Nm in the population.
    Nm_std (int, optional): The standard deviation of Nm in the population.
    Nf_avg (int, optional): The average Nf in the population.
    Nf_std (int, optional): The standard deviation of Nf in the population.
    covariance_avg (float, optional): The average covariance between loci.
    covariance_std (float, optional): The standard deviation of the covariance between loci.
    seed (int, optional): Random seed for reproducibility.

    Returns:
    dict: Dictionary with all possible genotypes and properties for 'male', 'female', and 'fitness'.

    Example dict (1 locus, 2 alleles):
    {
        (('A1', 'A1'),) {'male': 52, 'female': 46, 'fitness': 1.0, 'covariance': 0},
        (('A1', 'A2'),) {'male': 59, 'female': 82, 'fitness': 0.7, 'covariance': 0},
        (('A2', 'A2'),) {'male': 59, 'female': 44, 'fitness': 0.4, 'covariance': 0},
    }
    """

    if seed: random.seed(seed)

    loci_letters = [chr(ord('A') + i) for i in range(loci)]
    alleles_numbers = [f"{loci_letters[i]}{j+1}" for i in range(loci) for j in range(alleles)]

    def get_combinations(alleles):
        return list(itertools.combinations_with_replacement(alleles, 2))

    all_locus_combinations = [get_combinations(alleles_numbers[i*alleles:(i+1)*alleles]) for i in range(loci)]

    genotypes = list(itertools.product(*all_locus_combinations))

    genotype_data = defaultdict(lambda: {
        'Nm': int(random.gauss(Nm_avg, Nm_std)) if not Nm else Nm // len(genotypes), 
        'Nf': int(random.gauss(Nf_avg, Nf_std)) if not Nf else Nf // len(genotypes), 
        'fitness': round(random.uniform(0, 1), 1),
        'covariance': random.gauss(covariance_avg, covariance_std) if covariance_avg and covariance_std else 0
    })

    for genotype in genotypes:
        sorted_genotype = tuple(tuple(sorted(pair)) for pair in genotype)
        genotype_data[sorted_genotype]

    # Set fitness to 1 for a random genotype to emulate relative fitness
    random_genotype = random.choice(list(genotype_data.keys()))
    genotype_data[random_genotype]['fitness'] = 1.0

    return dict(genotype_data)

# # Example usage
# loci = 2
# alleles = 3
# Nm = 1000
# Nf = 1000
# seed = 42

# genotype_data = generate_genotype_data(loci, alleles, Nm, Nf, seed)

# for genotype, data in genotype_data.items():
#     print(genotype, data)