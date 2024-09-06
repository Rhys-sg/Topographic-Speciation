import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter


class TopographicMap:
    def __init__(self, width, height, name, smoothness=1.0, min_val=0, max_val=1, sum_val=None):
        self.width = width
        self.height = height
        self.name = name
        self.smoothness = smoothness
        self.min_val = min_val
        self.max_val = max_val
        self.sum_val = sum_val
        self.map_data = self.generate_random_topographical_map()

    def generate_random_topographical_map(self):
        """
        Generate a random topographical map with a normal distribution and normalize to [min_val, max_val].
        """
        map_data = np.random.normal(0.5, 0.2, (self.width, self.height))
        map_data = gaussian_filter(map_data, self.smoothness * ((self.width + self.height) // 2))
        map_data = (map_data - np.min(map_data)) / (np.max(map_data) - np.min(map_data))
        map_data = map_data * (self.max_val - self.min_val) + self.min_val
        
        # Adjust values to sum up to sum_val if provided
        if self.sum_val is not None:
            current_sum = np.sum(map_data)
            if current_sum != 0:
                scaling_factor = self.sum_val / current_sum
                map_data *= scaling_factor
        
        return map_data

    def __getitem__(self, idx):
        """
        Allows indexing using TopographicMap[x][y].
        """
        return self.map_data[idx]
    
    def __setitem__(self, idx, value):
        """
        Allows setting values using TopographicMap[x][y] = value.
        """
        self.map_data[idx] = value

    def plot_map(self):
        """
        Plot the random topographical map.
        """
        plt.figure(figsize=(10, 10))
        plt.imshow(self.map_data, cmap='gray', origin='upper')
        plt.colorbar()
        plt.title(f'{self.name} Topographical Map')
        plt.show()


# # Example Usage
# width = 100
# height = 100
# smoothness = 0.2
# min_val = 0
# max_val = 5
# sum_val = 10000

# topographical_map = TopographicMap(width, height, 'Example Map', smoothness, min_val, max_val, sum_val)
# topographical_map.plot_map()
