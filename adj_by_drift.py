import random

def adj_by_drift(curr_genotypes_data, max_drift, N, K):
    """
    Applies random genetic drift to the population, increasing or deceasing the genotype count up to max_drift.
    adj_drift is inversely relative to the population size with respect to carrying capacity.

    """

    adj_drift = max_drift * (1 - (N / K))
    for data in curr_genotypes_data.values():
        data['Nm'] = round(data['Nm'] * random.uniform(1-adj_drift, 1+adj_drift))
        data['Nf'] = round(data['Nf'] * random.uniform(1-adj_drift, 1+adj_drift))

    return curr_genotypes_data