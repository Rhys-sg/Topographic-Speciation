from calc_N import calc_N

def calc_avgerage_fitness(gens_genotype_data):
    """
    Calculate the average fitness of a population.

    """
    gens_avg_fitness = []
    for genotype_data in gens_genotype_data:
        avg_fitness = 0
        for data in genotype_data.values():
            avg_fitness += (data['Nm'] + data['Nf']) * data['fitness']

        N  = calc_N(genotype_data)

        gens_avg_fitness.append(avg_fitness / N)
        
    return gens_avg_fitness