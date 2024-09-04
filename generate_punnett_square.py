import itertools

def sorted_tuple(t):
    """
    Returns a sorted tuple of the input tuple.
    """
    return tuple(sorted(t))

def generate_punnett_square(genotype_A, genotype_B):
    """
    Generates a Punnett square from two genotypes with any number of loci and alleles.
    
    Args:
    - genotype_A, genotype_B: Tuples representing genotypes. Each tuple contains loci, where each locus is a list of alleles.
    
    Returns:
    - A dictionary representing the Punnett square, where keys are sorted combinations of alleles and values are counts of the corresponding genotypes.
    """
    
    # Prepare the loci for each genotype
    loci_combinations = []
    for locus1, locus2 in zip(genotype_A, genotype_B):
        # Create all possible combinations for this locus
        combinations = [sorted_tuple((allele1, allele2)) for allele1, allele2 in itertools.product(locus1, locus2)]
        loci_combinations.append(combinations)
    
    # Generate all possible genotype combinations
    punnett_square = {}
    for combination in itertools.product(*loci_combinations):
        # Sort alleles in each locus to handle equivalent combinations
        sorted_combination = tuple(sorted_tuple(comb) for comb in combination)
        if sorted_combination not in punnett_square:
            punnett_square[sorted_combination] = 0
        punnett_square[sorted_combination] += 1

    # normalize the punnett square
    total = sum(punnett_square.values())
    for key in punnett_square:
        punnett_square[key] /= total

    return punnett_square