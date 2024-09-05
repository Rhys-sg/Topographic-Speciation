import itertools
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from topographical_map import TopographicMap
from cell import Cell

# popgen functions
from calc_N import calc_N
from calc_next_N import calc_next_N
from calc_N_sub import calc_N_sub

from adj_by_fitness import adj_by_fitness
from adj_by_drift import adj_by_drift
from adj_by_mutation import adj_by_mutation

from calc_next_genotypes_data import calc_next_genotypes_data

from animator import Animator

class TopographicalSpeciation:
    """
    A class to simulate the speciation of subpopulations in the cells of a 2d map.
    The carrying capcity and the fitness of each genotype is determined by a randomly generated topographical map.

    TODO:
    - Implement gene flow
    - Implement covariance
    - Implement phisical barriers

    """
    def __init__(self, width, height, N, K, smoothness=1.0, loci=1, alleles=2, growth_rate=0, max_drift=0, mutation_rate=None):
        self.width = width
        self.height = height
        self.N = N
        self.K = K
        self.smoothness = smoothness
        self.loci = loci
        self.alleles = alleles

        self.growth_rate = growth_rate
        self.max_drift = max_drift
        self.mutation_rate = mutation_rate

        self.generations = 0

        self.genotypes = self.generate_genotypes(self.alleles, self.loci)
        self.cells = self.generate_cells()
        self.topographical_fitnesses = self.generate_topographical_fitnesses()

        self.animator = Animator()

    def generate_cells(self):
        cells = [[None for _ in range(self.width)] for _ in range(self.height)]

        carrying_capacity_map = TopographicMap(self.width, self.height, 'Carrying Capacity', self.smoothness)

        init_cell_N = self.N / (self.width * self.height)
        total = sum([carrying_capacity_map.map_data[y][x] for x in range(self.width) for y in range(self.height)])
        for x in range(self.width):
            for y in range(self.height):
                carrying_capacity = (carrying_capacity_map.map_data[y][x] / total) * self.K
                cells[y][x] = Cell(x, y, init_cell_N/2, init_cell_N/2, self.genotypes, carrying_capacity)

        return cells
    
    def generate_topographical_fitnesses(self):
        topographical_fitnesses = {}
        for genotype in self.genotypes:
            topographical_fitnesses[genotype] = TopographicMap(self.width, self.height, genotype, self.smoothness)
        return topographical_fitnesses
    
    def generate_genotypes(self, alleles, loci):

        loci_letters = [chr(ord('A') + i) for i in range(loci)]
        alleles_numbers = [f"{loci_letters[i]}{j+1}" for i in range(loci) for j in range(alleles)]

        def get_combinations(alleles):
            return list(itertools.combinations_with_replacement(alleles, 2))

        all_locus_combinations = [get_combinations(alleles_numbers[i*alleles:(i+1)*alleles]) for i in range(loci)]

        return list(itertools.product(*all_locus_combinations))
    
    def run(self, generations, bottleneck_yr=None, bottleneck_N=None):
        self.generations += generations

        for i in range(generations):
            self.calc_generation(i, bottleneck_yr, bottleneck_N)
        
        return self.cells

    def calc_generation(self, i, bottleneck_yr=None, bottleneck_N=None):
        for x in range(self.width):
            for y in range(self.height):
                cell = self.cells[y][x]

                curr_genotypes_data = cell.gens_genotype_data[-1]

                next_N = calc_next_N(calc_N(curr_genotypes_data), self.growth_rate, cell.carrying_capacity)

                # Apply evolutionary forces to genotypes in the current generation
                curr_genotypes_data = adj_by_fitness(curr_genotypes_data, self.topographical_fitnesses, x, y)
                curr_genotypes_data = adj_by_drift(curr_genotypes_data, self.max_drift, calc_N(curr_genotypes_data), cell.carrying_capacity)

                # If bottleneck, adjust population size, apply drift
                if bottleneck_yr == i and next_N > bottleneck_N:
                    curr_genotypes_data = adj_by_drift(curr_genotypes_data, 1-bottleneck_N/next_N, bottleneck_N, cell.carrying_capacity)
                    next_N = bottleneck_N

                # Calculate the next generation
                curr_genotypes_data = calc_next_genotypes_data(curr_genotypes_data, next_N)
                curr_genotypes_data = adj_by_mutation(curr_genotypes_data, self.mutation_rate)

                cell.gens_genotype_data.append(curr_genotypes_data)


    def animate_gens_genotypes(self, genotypes=None):
        if not genotypes:
            genotypes = self.genotypes

        self.animator.animate_gens_genotypes(genotypes, self.cells, self.width, self.height, self.generations)
    

    def animate_gens_genotype_prevalence(self, genotypes=None, gradient=True):
        if not genotypes:
            genotypes = self.genotypes
        
        self.animator.animate_gens_genotype_prevalence(genotypes, self.cells, self.width, self.height, self.generations, gradient)



if __name__ == "__main__":
    width = 100
    height = 100
    N = 20 * width * height
    K = 200 * width * height
    smoothness = 0.1
    loci = 1
    alleles = 2

    growth_rate = 1.1
    mutation_rate=0.01
    # max_drift=0.1
    max_drift=0

    ts = TopographicalSpeciation(width, height, N, K, smoothness, loci, alleles, growth_rate, max_drift, mutation_rate)

    ts.run(20)
    ts.animate_gens_genotypes()
    ts.animate_gens_genotype_prevalence()
