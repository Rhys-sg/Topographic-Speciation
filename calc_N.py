def calc_N(curr_genotypes_data):
    """
    Counts the poopulation size for males and females
    
    """
    return sum([data['Nm'] + data['Nf'] for data in curr_genotypes_data.values()])