import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the tab-separated txt file (adjust path as needed)
df = pd.read_csv("/Users/milicaurosevic/Desktop/TP53_test/Benchmark/summary_counters_mousee.txt", sep="\t")

genes = df['Gene']

# Columns with percentage counts
counts = ['Count1(%)', 'Count2(%)', 'Count3(%)', 'Count4(%)']

data = df[counts].values

x = np.arange(len(genes))
width = 0.7

# Colors from light to dark for increasing number of supporting tools
colors = ['#a6cee3', '#1f78b4', '#b2df8a', '#33a02c']

fig, ax = plt.subplots(figsize=(14,7))

bottom = np.zeros(len(genes))

# Plot stacked bars
for i in range(len(counts)):
    ax.bar(x, data[:, i], width, bottom=bottom, color=colors[i], 
           label=f'Supported by at least {i+1} tool{"s" if i > 0 else ""}')
    bottom += data[:, i]

ax.set_xticks(x)
ax.set_xticklabels(genes, rotation=45, ha='right')

ax.set_ylabel('Percentage of new paths discovered (%)')
ax.set_xlabel('Gene')
ax.set_title('The proportion of new paths discovered - mouse')

ax.set_ylim(0, 100)  # Set y-axis scale from 0 to 100%

ax.legend(title='Support Level', bbox_to_anchor=(1.05, 1), loc='upper left')

plt.tight_layout()
plt.savefig("stacked_barplot_mouse.png", dpi=300)
plt.show()
