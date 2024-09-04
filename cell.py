import random
from collections import defaultdict

class Cell:
    def __init__(self, x, y, Nm, Nf, genotypes, carrying_capacity):
        self.x = x
        self.y = y
        self.Nm = Nm
        self.Nf = Nf
        self.genotypes = genotypes
        self.carrying_capacity = carrying_capacity

        self.gens_genotype_data = [self.generate_genotype_data()]
        
    def generate_genotype_data(self, covariance_avg=None, covariance_std=None):

        genotype_data = defaultdict(lambda: {
            'Nm':  self.Nm // len(self.genotypes), 
            'Nf':  self.Nf // len(self.genotypes), 
            'covariance': random.gauss(covariance_avg, covariance_std) if covariance_avg and covariance_std else 0
        })

        for genotype in self.genotypes:
            sorted_genotype = tuple(tuple(sorted(pair)) for pair in genotype)
            genotype_data[sorted_genotype]

        return dict(genotype_data)
    
