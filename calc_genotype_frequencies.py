from calc_N import calc_N

def calc_genotype_frequencies(gens_genotype_data):
    """
    Calculate genotype frequencies from genotype counts.
    """
    return [{genotype: (data['Nm'] + data['Nf']) / calc_N(genotype_data) for genotype, data in genotype_data.items()} for genotype_data in gens_genotype_data]
