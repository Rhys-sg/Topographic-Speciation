import random

def adj_by_mutation(genotype_data, mutation_rate):
    mutated_data = {k: v.copy() for k, v in genotype_data.items()}

    for genotype_a, data_a in genotype_data.items():
        for genotype_b, data_b in genotype_data.items():
            if genotype_a == genotype_b:
                continue
            distance = calc_distance(genotype_a, genotype_b)
            mutation_prob = mutation_rate ** distance

            total_individuals = data_a['Nm'] + data_a['Nf']
            if total_individuals > 0:
                N_mutations = int((total_individuals * mutation_prob) / (len(genotype_data) * 2))
                
                # Ensure mutations do not result in negative counts
                N_mutations = min(N_mutations, data_a['Nm'])
                N_mutations = min(N_mutations, data_a['Nf'])

                mutated_data[genotype_a]['Nm'] -= N_mutations
                mutated_data[genotype_a]['Nf'] -= N_mutations
                mutated_data[genotype_b]['Nm'] += N_mutations
                mutated_data[genotype_b]['Nf'] += N_mutations

    return mutated_data

def calc_distance(genotype_a, genotype_b):
    distance = 0
    for locus in range(len(genotype_a)):
        for allele in range(len(genotype_a[locus])):
            if genotype_a[locus][allele] != genotype_b[locus][allele]:
                distance += 1
    return distance
