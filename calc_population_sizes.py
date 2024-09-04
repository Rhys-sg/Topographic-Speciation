from calc_N import calc_N

def calc_population_sizes(gens_genotype_data):
    
    return [calc_N(genotype_data) for genotype_data in gens_genotype_data]
