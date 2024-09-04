def calc_allele_counts(gens_genotype_data):
    """
    Calculate allele counts for multiple loci with given the male/female genotype counts.

    """
    gens_allele_counts = []
    for genotype_data in gens_genotype_data:
        allele_counts = {}

        for genotype, data in genotype_data.items():
            for i, locus in enumerate(genotype):
                for allele in locus:
                    if allele not in allele_counts:
                        allele_counts[allele] = 0
                    allele_counts[allele] += data['Nm'] + data['Nf']
                    
        gens_allele_counts.append(allele_counts)

    return gens_allele_counts