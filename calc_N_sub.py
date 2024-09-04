def calc_N_sub(gens_genotype_data, sub_pop):
    """
    Counts the poopulation size for males or females
    
    """
    return [sum([data[sub_pop] for data in genotype_data.values()]) for genotype_data in gens_genotype_data]
