def append_gene_flow_cells(cells, topographical_gene_flow_distance):
    def update_gene_flow(x, y, Nm_or_Nf, new_N, rate, max_distance):
        if rate > 0 and Nm_or_Nf > 0 and max_distance > 0:
            for i in range(max(0, x-max_distance), min(len(cells), x+max_distance+1)):
                for j in range(max(0, y-max_distance), min(len(cells[0]), y+max_distance+1)):
                    distance = abs(i-x) + abs(j-y)
                    new_N[i][j] += (rate * Nm_or_Nf) / distance if distance > 0 else Nm_or_Nf * (1-rate)
        else:
            new_N[x][y] += Nm_or_Nf

    new_Nm = {genotype: [[0] * len(cells[0]) for _ in range(len(cells))] for genotype in cells[0][0].genotypes}
    new_Nf = {genotype: [[0] * len(cells[0]) for _ in range(len(cells))] for genotype in cells[0][0].genotypes}

    for x in range(len(cells)):
        for y in range(len(cells[0])):
            N = sum(cells[x][y].gens_genotype_data[-1][genotype]['Nm'] + cells[x][y].gens_genotype_data[-1][genotype]['Nf'] for genotype in cells[x][y].genotypes)
            K = cells[x][y].carrying_capacity

            if N > K and N > 0:
                rate = (N - K) / N
            else:
                rate = 0

            for genotype in cells[x][y].genotypes:
                max_distance = int(topographical_gene_flow_distance[genotype].map_data[y][x])
                Nm = cells[x][y].gens_genotype_data[-1][genotype]['Nm']
                Nf = cells[x][y].gens_genotype_data[-1][genotype]['Nf']

                update_gene_flow(x, y, Nm, new_Nm[genotype], rate, max_distance)
                update_gene_flow(x, y, Nf, new_Nf[genotype], rate, max_distance)

    for x in range(len(cells)):
        for y in range(len(cells[0])):
            genotype_data = {
                genotype: {
                    'Nm': new_Nm[genotype][x][y],
                    'Nf': new_Nf[genotype][x][y],
                    'covariance': cells[x][y].gens_genotype_data[0][genotype]['covariance'],
                }
                for genotype in cells[x][y].genotypes
            }
            cells[x][y].gens_genotype_data.append(genotype_data)

    return cells



def append_cells(cells):
    for x in range(len(cells)):
        for y in range(len(cells[0])):
            cells[x][y].gens_genotype_data.append(cells[x][y].gens_genotype_data[-1])
    return cells