import itertools
from collections import defaultdict
import random

from topographical_map import TopographicMap


class TopographicalSpeciation:
    def __init__(self, width, height, N, smoothness=1.0, loci=1, alleles=2):
        self.width = width
        self.height = height
        self.N = N
        self.smoothness = smoothness
        self.loci = loci
        self.alleles = alleles

        self.genotypes = self.generate_genotypes(self.alleles, self.loci)
        self.carrying_capacity = TopographicMap(self.width, self.height, 'Carrying Capacity', self.smoothness)
        self.cells = self.generate_cells()
        self.topographical_fitnesses = self.generate_topographical_fitnesses()

    def generate_cells(self):
        cells = [[None for _ in range(self.width)] for _ in range(self.height)]

        init_cell_N = self.N / (self.width * self.height)
        for x in range(self.width):
            for y in range(self.height):
                cells[y][x] = Cell(x, y, init_cell_N/2, init_cell_N/2, self.genotypes, self.carrying_capacity.map_data[y][x])

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
        for i in range(generations):
            self.calc_generation(i, bottleneck_yr, bottleneck_N)
        
        return self.cells

    def calc_generation(self, i, bottleneck_yr=None, bottleneck_N=None):
        pass


class Cell:
    def __init__(self, x, y, Nm, Nf, genotypes, carrying_capacity):
        self.x = x
        self.y = y
        self.Nm = Nm
        self.Nf = Nf
        self.genotypes = genotypes
        self.carrying_capacity = carrying_capacity

        self.genotype_data = self.generate_genotype_data()

    def generate_genotype_data(self, covariance_avg=None, covariance_std=None):

        genotype_data = defaultdict(lambda: {
            'Nm':  self.Nm // len(self.genotypes), 
            'Nf':  self.Nf // len(self.genotypes), 
            'fitness': round(random.uniform(0, 1), 1),
            'covariance': random.gauss(covariance_avg, covariance_std) if covariance_avg and covariance_std else 0
        })

        for genotype in self.genotypes:
            sorted_genotype = tuple(tuple(sorted(pair)) for pair in genotype)
            genotype_data[sorted_genotype]

        # Set fitness to 1 for a random genotype to emulate relative fitness
        random_genotype = random.choice(list(genotype_data.keys()))
        genotype_data[random_genotype]['fitness'] = 1.0

        return dict(genotype_data)


if __name__ == "__main__":
    width = 50
    height = 50
    N = 100
    smoothness = 3
    loci = 1
    alleles = 2

    topographical_speciation = TopographicalSpeciation(width, height, N, smoothness, loci, alleles)
    # for each in topographical_speciation.topographical_fitnesses:
    #     each.plot_map()

    # print(topographical_speciation.cells[0][0].genotype_data)
    # for cell in topographical_speciation.cells:
    #     for genotype, data in cell.items():
    #         print(f"{genotype}: {data}")

    topographical_speciation.carrying_capacity.plot_map()