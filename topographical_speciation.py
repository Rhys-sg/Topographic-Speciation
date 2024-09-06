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

from append_cells import append_cells, append_gene_flow_cells

from animator import Animator

class TopographicalSpeciation:
    """
    A class to simulate the speciation of subpopulations in the cells of a 2d map.
    The carrying capcity and the fitness of each genotype is determined by a randomly generated topographical map.

    TODO:
    - Implement covariance
    - Implement physical barriers
        - Use a topographical map to represent how far a population can move from that cell

    """
    def __init__(self):
        self.width = None
        self.height = None
        self.N = None
        self.K = None
        self.min_cell_K = 0

        self.smoothness = 1.0
        self.loci = 1
        self.alleles = 2

        self.growth_rate = 0
        self.max_drift = 0
        self.max_gene_flow_distance = 0
        self.mutation_rate = None

        self.generations = 0

        self.animator = Animator()

    def generate_fields(self):
        if self.width is None or self.height is None or self.N is None or self.K is None:
            raise ValueError("Width, height, N, and K must be set before generating fields.")
        
        self.topographical_carrying_capacity = TopographicMap(self.width,
                                                              self.height,
                                                              'Carrying Capacity',
                                                              self.smoothness,
                                                              min_val=self.min_cell_K,
                                                              sum_val=self.K)
        self.genotypes = self.generate_genotypes(self.alleles, self.loci)
        self.cells = self.generate_cells()

        self.topographical_fitnesses = self.generate_topographical_map_for_genotype()
        self.topographical_gene_flow_distance = self.generate_topographical_map_for_genotype(max_val=self.max_gene_flow_distance)

    def generate_cells(self):
        cells = [[None for _ in range(self.width)] for _ in range(self.height)]

        cell_N = self.N / (self.width * self.height)
        for x in range(self.width):
            for y in range(self.height):
                carrying_capacity = self.topographical_carrying_capacity[y][x]
                cells[y][x] = Cell(x, y, cell_N/2, cell_N/2, self.genotypes, carrying_capacity)

        return cells
    
    def generate_topographical_map_for_genotype(self, min_val=0, max_val=1, sum_val=None):
        topographical_maps = {}
        for genotype in self.genotypes:
            topographical_maps[genotype] = TopographicMap(self.width,
                                                          self.height,
                                                          genotype,
                                                          self.smoothness,
                                                          min_val=min_val,
                                                          max_val=max_val,
                                                          sum_val=sum_val)
        return topographical_maps
    
    def generate_genotypes(self, alleles, loci):

        loci_letters = [chr(ord('A') + i) for i in range(loci)]
        alleles_numbers = [f"{loci_letters[i]}{j+1}" for i in range(loci) for j in range(alleles)]

        def get_combinations(alleles):
            return list(itertools.combinations_with_replacement(alleles, 2))

        all_locus_combinations = [get_combinations(alleles_numbers[i*alleles:(i+1)*alleles]) for i in range(loci)]

        return list(itertools.product(*all_locus_combinations))
    
    def run(self, generations, bottleneck_yr=None, bottleneck_N=None):
        self.generate_fields()

        self.generations += generations

        for i in range(generations):
            print(f"Generation {i}")
            self.calc_generation(i, bottleneck_yr, bottleneck_N)
        
        return self.cells

    def calc_generation(self, i, bottleneck_yr=None, bottleneck_N=None):

        # Initialize the next generation, apply gene flow
        if self.max_gene_flow_distance > 0:
            self.cells = append_gene_flow_cells(self.cells, self.topographical_gene_flow_distance)
        else:
            self.cells = append_cells(self.cells) 

        # For each cell in the map
        for x in range(self.width):
            for y in range(self.height):
                cell = self.cells[y][x]

                curr_genotypes_data = cell.gens_genotype_data[-1]

                # next_N = min(calc_N(curr_genotypes_data), cell.carrying_capacity)
                next_N = calc_next_N(calc_N(curr_genotypes_data), self.growth_rate, cell.carrying_capacity)


                # Apply evolutionary forces to genotypes in the current generation
                curr_genotypes_data = adj_by_fitness(curr_genotypes_data, self.topographical_fitnesses, x, y)
                curr_genotypes_data = adj_by_drift(curr_genotypes_data, self.max_drift, calc_N(curr_genotypes_data), cell.carrying_capacity)

                # If bottleneck, adjust population size, apply drift
                if bottleneck_yr == i and next_N > bottleneck_N:
                    curr_genotypes_data = adj_by_drift(curr_genotypes_data, 1-bottleneck_N/next_N, bottleneck_N, cell.carrying_capacity)
                    next_N = bottleneck_N

                # Calculate the next generation
                curr_genotypes_data = calc_next_genotypes_data(curr_genotypes_data, next_N * self.growth_rate)
                curr_genotypes_data = adj_by_mutation(curr_genotypes_data, self.mutation_rate)

                # Update the current cell's last genotypes data
                cell.gens_genotype_data[-1] = curr_genotypes_data


    def animate_gens_genotypes(self, genotypes=None):
        if not genotypes:
            genotypes = self.genotypes

        self.animator.animate_gens_genotypes(genotypes, self.cells, self.width, self.height, self.generations)
    

    def animate_gens_genotype_prevalence(self, genotypes=None, gradient=True):
        if not genotypes:
            genotypes = self.genotypes
        
        self.animator.animate_gens_genotype_prevalence(genotypes, self.cells, self.width, self.height, self.generations, gradient)



if __name__ == "__main__":
    ts = TopographicalSpeciation()

    ts.width = 20
    ts.height = 20
    ts.N = 20 * ts.width * ts.height
    ts.K = 200 * ts.width * ts.height
    ts.min_cell_K = 10
    ts.smoothness = 0.15
    ts.loci = 1
    ts.alleles = 2
    ts.growth_rate = 1.1
    ts.max_drift=0
    ts.max_gene_flow_distance=3
    ts.mutation_rate=0.01

    ts.run(30)
    ts.animate_gens_genotypes()
    # ts.animate_gens_genotype_prevalence()
