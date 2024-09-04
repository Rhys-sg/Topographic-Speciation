from calc_allele_counts import calc_allele_counts

def calc_allele_frequencies(gens_genotype_data):

    gens_allele_counts = calc_allele_counts(gens_genotype_data)
    gens_allele_freq = []

    for allele_counts in gens_allele_counts:
        N = sum(allele_counts.values())
        locus_freqs = {allele: count / N for allele, count in allele_counts.items()}
        gens_allele_freq.append(locus_freqs)
    

    return gens_allele_freq