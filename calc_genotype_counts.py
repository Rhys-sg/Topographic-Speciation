def calc_genotype_counts(gens_genotype_data):
    """
    Calculate genotype counts from genotype data.
    
    """

    return [{genotype : data['Nm'] + data['Nf'] for genotype, data in genotype_data.items()} for genotype_data in gens_genotype_data]