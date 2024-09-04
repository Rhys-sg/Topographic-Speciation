def adj_by_fitness(curr_genotypes_data, topographical_fitnesses, x, y):
    """
    Adjust genotype counts by their fitness.

    """
    for genotype, data in curr_genotypes_data.items():
        data['Nm'] = round(data['Nm'] * topographical_fitnesses[genotype][x][y])
        data['Nf'] = round(data['Nf'] * topographical_fitnesses[genotype][x][y])
    
    return curr_genotypes_data