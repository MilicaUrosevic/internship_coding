import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Учитавање CSV фајла (путања прилагодити по потреби)
df = pd.read_csv("/Users/milicaurosevic/Desktop/TP53_test/Benchmark/results_final.csv", sep=",")

genes = df['Gene']

cols_first_bar = [
    'with_4_suffix_and_gene', 
    'with_3_suffix_and_gene', 
    'with_2_suffix_and_gene', 
    'with_1_suffix_and_gene'
]

cols_second_bar = [
    'without_4_suffix_and_gene',
    'without_3_suffix_and_gene',
    'without_2_suffix_and_gene',
    'without_1_suffix_and_gene'
]

legend_labels_first = [
    's-exons supported by at least 1 tool and the annotation', 
    's-exons supported by at least 2 tools and the annotation', 
    's-exons supported by at least 3 tools and the annotation', 
    's-exons supported by at least 4 tools and the annotation'
]

legend_labels_second = [
    's-exons gained by at least 1 tool',
    's-exons gained by at least 2 tools',
    's-exons gained by at least 3 tools',
    's-exons gained by at least 4 tools'
]

# Боје — плаве нијансе за прву групу (од најтамније ка светлијој)
colors_first = ['#084594', '#2171b5', '#4292c6', '#6baed6']

# Боје — наранџасте нијансе за другу групу (од најтамније ка светлијој)
colors_second = ['#a63603', '#e6550d', '#fd8d3c', '#fdae6b']

# Шире размаци између гена (x оса)
spacing = 1.5  
x = np.arange(len(genes)) * spacing  

# Ширина барова сада шира
width = 0.45  

fig, ax = plt.subplots(figsize=(14, 7))

def plot_stacked_bar(ax, x_positions, data, offset, colors):
    bottoms = np.zeros(len(data[0]))
    bars = []
    for i, col in enumerate(data):
        bar = ax.bar(x_positions + offset, col, width, bottom=bottoms, color=colors[i])
        bottoms += col
        bars.append(bar)
    return bars

data_first = [df[col].values for col in cols_first_bar]
data_second = [df[col].values for col in cols_second_bar]

bars_first = plot_stacked_bar(ax, x, data_first, -width/2, colors_first)
bars_second = plot_stacked_bar(ax, x, data_second, width/2, colors_second)

ax.set_xticks(x)
ax.set_xticklabels(genes, rotation=45, ha='right')

ax.set_ylabel('Percentage of s-exon representation')
ax.set_xlabel('Gene ID')
ax.set_title('Stacked Barplot of s-exon representaton per gene, mouse')

ax.set_ylim(0, 1.0)

legend_labels_first_rev = legend_labels_first[::-1]
legend_labels_second_rev = legend_labels_second[::-1]

legend_bars = [bars_first[i][0] for i in range(4)] + [bars_second[i][0] for i in range(4)]
legend_labels = legend_labels_first_rev + legend_labels_second_rev

ax.legend(legend_bars, legend_labels, title="Percentage of supported (left) and gained (right) s-exons ", bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("barplot_s_mouse.png", dpi=300)
plt.show()
