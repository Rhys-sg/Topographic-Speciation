import plotly.graph_objects as go
import plotly.colors as pc

def format_genotype_key(genotype):
    return ', '.join(''.join(pair) for pair in genotype)

def create_plot(data_type, gens_data, y_label, mode='lines'):
    color_cycle = pc.qualitative.Plotly

    # Collect all keys for color assignment
    all_keys = set()
    for gen in gens_data:
        if isinstance(gen, dict):
            all_keys.update(gen.keys())
        else:
            all_keys.add(data_type)

    colors = {key: color_cycle[i % len(color_cycle)] for i, key in enumerate(all_keys)}

    # Process data
    def collect_data(gens_data):
        data = {}
        for gens_index, gen_data in enumerate(gens_data):
            if isinstance(gen_data, dict):
                for key, value in gen_data.items():
                    if key not in data:
                        data[key] = [0] * len(gens_data)
                    data[key][gens_index] = value
            else:
                if data_type not in data:
                    data[data_type] = [0] * len(gens_data)
                data[data_type][gens_index] = gen_data
        return data

    # Collect data
    data = collect_data(gens_data)

    # Create traces
    traces = []
    for key, values in data.items():
        trace = go.Scatter(
            x=list(range(len(values))),
            y=values,
            mode=mode,
            name=format_genotype_key(key) if isinstance(key, tuple) else key,
            line=dict(color=colors.get(key, 'black')),
            hoverinfo='text',
            text=[f'{format_genotype_key(key) if isinstance(key, tuple) else key}: {value}' for value in values]
        )
        traces.append(trace)

    # Create figure
    fig = go.Figure(data=traces)

    # Update layout
    fig.update_layout(
        title=f'{data_type.replace("_", " ").title()} Over Generations',
        xaxis_title='Generation',
        yaxis_title=y_label,
        hovermode='closest',
        height=600
    )

    # Update y-axis range for frequency plots
    if 'Frequency' in y_label:
        fig.update_yaxes(range=[0, 1])

    return fig
