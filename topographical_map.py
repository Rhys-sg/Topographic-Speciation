import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter

class TopographicMap:
    def __init__(self, width, height, genotype, smoothness=1.0):
        self.width = width
        self.height = height
        self.genotype = genotype
        self.smoothness = smoothness
        self.map_data = self.generate_random_topographical_map()

    def generate_random_topographical_map(self):
        """
        Generate a random topographical map for a fitness with a normal distribution from 0 to 1.
        """
        map_data = np.random.normal(0.5, 0.2, (self.width, self.height))
        map_data = gaussian_filter(map_data, self.smoothness)
        map_data = (map_data - np.min(map_data)) / (np.max(map_data) - np.min(map_data))
        
        return map_data

    def plot_map(self):
        """
        Plot the random topographical map.
        """
        plt.figure(figsize=(10, 10))
        plt.imshow(self.map_data, cmap='gray', origin='upper')
        plt.colorbar()
        plt.title(f'{self.genotype} Topographical Map')
        plt.show()


# # Parameters
# width = 512
# height = 512
# smoothness = 50

# # Create a topographical map
# topographical_map = topographical_fitness(width, height, smoothness)
# topographical_map.plot_map()