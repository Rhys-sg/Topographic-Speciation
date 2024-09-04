def calc_next_N(N, r, K):
    """
    Calculates the size of the next population based on Verhulst's model:
    - dN/dt = rate of population change
    -     r = maximum population growth rate
    -     N = current population size
    -     K = population carrying capacity
    """

    if not K:
        dN = N * r
    else:
        dN = r * N * (1 - (N / K))

    # Ensure the population size never goes negative
    next_N = max(1, N + dN)

    return next_N
