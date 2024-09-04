import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

class Animator:
    def animate_gens_genotypes(self, genotypes, cells, width, height, generations):

        num_genotypes = len(genotypes)
        fig, axs = plt.subplots(1, num_genotypes, figsize=(10 * num_genotypes, 10), sharex=True, sharey=True)
        
        if num_genotypes == 1:
            axs = [axs]

        # Placeholder for the initial image
        im = None

        def update(i):
            nonlocal im
            for idx, genotype in enumerate(genotypes):
                map_data = [[None for _ in range(width)] for _ in range(height)]
                for x in range(width):
                    for y in range(height):
                        cell = cells[y][x]
                        N = cell.gens_genotype_data[i][genotype]['Nm'] + cell.gens_genotype_data[i][genotype]['Nf']
                        percentage = N / cell.carrying_capacity if cell.carrying_capacity != 0 else 0
                        map_data[y][x] = percentage

                axs[idx].clear()
                im = axs[idx].imshow(map_data, cmap='gray', origin='upper', vmin=0, vmax=1)
                axs[idx].set_title(f'{genotype} at Generation {i}')

            return [im] * num_genotypes

        ani = animation.FuncAnimation(fig, update, frames=generations, interval=300, repeat=True)

        dummy_plot = axs[0].imshow([[0]], cmap='gray', origin='upper', vmin=0, vmax=1)
        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        cbar = plt.colorbar(dummy_plot, cax=cbar_ax)
        cbar.set_label('Percentage')

        plt.show()
    

    def animate_gens_genotype_prevalence(self, genotypes, cells, width, height, generations, gradient):

        fig, ax = plt.subplots(figsize=(10, 10))
        num_genotypes = len(genotypes)
        colors = plt.cm.get_cmap('tab20', num_genotypes)

        def update(i):
            map_data = np.zeros((height, width), dtype=int)
            for x in range(width):
                for y in range(height):
                    cell = cells[y][x]
                    max_percentage = -1
                    most_prevalent_genotype = None

                    for genotype in genotypes:
                        genotype_data = cell.gens_genotype_data[i].get(genotype, {'Nm': 0, 'Nf': 0})
                        N = genotype_data['Nm'] + genotype_data['Nf']
                        percentage = N / cell.carrying_capacity if cell.carrying_capacity != 0 else 0

                        if percentage > max_percentage:
                            max_percentage = percentage
                            most_prevalent_genotype = genotype
                    
                    # Set the map data to the index of the most prevalent genotype
                    if most_prevalent_genotype is not None:
                        map_data[y, x] = genotypes.index(most_prevalent_genotype)
            
            ax.clear()
            cax = ax.imshow(map_data, cmap=colors, origin='upper', vmin=0, vmax=num_genotypes-1)
            ax.set_title(f'Most Prevalent Genotype at Generation {i}')
            return cax,

        ani = animation.FuncAnimation(fig, update, frames=generations, interval=300, repeat=True)
        
        cbar = plt.colorbar(ax.imshow(np.zeros((1, 1)), cmap=colors, origin='upper', vmin=0, vmax=num_genotypes-1), ax=ax)
        cbar.set_ticks(np.arange(num_genotypes))
        cbar.set_ticklabels(genotypes)
        cbar.set_label('Genotype')

        plt.show()