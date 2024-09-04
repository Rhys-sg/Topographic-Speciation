# Wikipedia: https://en.wikipedia.org/wiki/Effective_population_size

def calc_Ne_sex_ratios(Nm, Nf):
    """
    Calculate effective population size (Ne) based on Non-Fisherian sex-ratios
    
    Args:
    - Nm: Number of males
    - Nf: Number of females
    
    Returns:
    - Ne: Effective population size
    """
    Ne = (4 * Nm * Nf) / (Nm + Nf)
    return Ne


def calc_Ne_size_variation(N):
    """
    Calculate effective population size (Ne) based on variation for all previous generations.
    
    Args:
    - N: List of population sizes for all previous generations

    Returns:
    - harmonic_mean: Harmonic mean of population sizes
    """
    reciprocal_sum = sum(1.0 / N_t for N_t in N)

    return len(N) / reciprocal_sum

# Hutch equations

def calc_Ne_change_in_allele_freq(allele_freq, N, t):
    """
    Calculate effective population size (Ne) based on the statistical measure of
    allelic variation (sigma) at a given generation.
    
    If N is really large, then the part in the parentheses becomes 1-0 and the part
    in the brackets becomes 1-1, or zero, making the equation to zero. When this is true,
    there is no variance and no bouncing around of allele frequencies. 

    Args:
    - allele_freq: List of allele frequencies
    - N: Population size
    - t: Generation number

    Returns:
    - Ne: Effective population size based on allele variation

    """

    sigma = sum(p * (1 - p) for p in allele_freq) * ((1-(1/(2*N))) ** t)
    if sigma == 0:
        return N
    return N * (1-sigma)


def calc_Ne_inbreeding(N, t):
    """
    Calculate effective population size (Ne) based on inbreeding coefficient (F).
    This is according to Hutch, there is another equation on Wikipedia defining it
    by how the average inbreeding coefficient changes from one generation to the next

    As N increases, the part in the parentheses reduces to 1-0 and the equation approaches zero. When
    this is true, there is no increase in identity by descent due to drift.

    Args:
    - N: Population size
    - t: Generation number
    """

    F = 1 - ((1 - (1 / N)) ** t)
    return N * (1 - F)



def calc_Ne_over_generations(gens_Nm, gens_Nf, gens_allele_freqs):
    N = [Nm + Nf for Nm, Nf in zip(gens_Nm, gens_Nf)]
    gens_allele_freq_list = [list(freq_dict.values()) for freq_dict in gens_allele_freqs]

    gens_Ne = []
    for i in range(len(N)):
        gen_Ne = {
            'Ne (Sex Ratios)': calc_Ne_sex_ratios(gens_Nm[i], gens_Nf[i]),
            'Ne (Size Variation)': calc_Ne_size_variation(N[:i+1]),
            'Ne (Allele Variation)': calc_Ne_change_in_allele_freq(gens_allele_freq_list[i], N[i], i + 1),
            'Ne (Inbreeding)': calc_Ne_inbreeding(N[i], i + 1)
        }
        gens_Ne.append(gen_Ne)

    return gens_Ne
