import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df_both = pd.read_csv("support_events_both_mouse.csv")
df_alt = pd.read_csv("support_events_alt_mouse.csv")

for df in [df_both, df_alt]:
    df['1'] = df['1']
    df['2'] = df['1'] + df['2']
    df['3'] = df['2'] + df['3']
    df['4'] = df['3'] + df['4']

genes = df_both['gene']
x = np.arange(len(genes))
width = 0.4

colors_both = ['#c6dbef', '#9ecae1', '#6baed6', '#3182bd']
colors_alt  = ['#fdd0a2', '#fdae6b', '#fd8d3c', '#e6550d']

fig, ax = plt.subplots(figsize=(14, 7))

bottom = np.zeros(len(genes))
for i, col in enumerate(['1', '2', '3', '4']):
    heights = df_both[col] - bottom
    label = f'{i+1} tool{"s" if i>0 else ""} canonical + alternative'
    ax.bar(x - width/2, heights, width, bottom=bottom, color=colors_both[i], label=label)
    bottom = df_both[col]

bottom = np.zeros(len(genes))
for i, col in enumerate(['1', '2', '3', '4']):
    heights = df_alt[col] - bottom
    label = f'{i+1} tool{"s" if i>0 else ""} alternative'
    ax.bar(x + width/2, heights, width, bottom=bottom, color=colors_alt[i], label=label)
    bottom = df_alt[col]

ax.set_ylabel('Percentage of supported events per tool')
ax.set_title('Supported splicing events in genes, Mus musculus')
ax.set_xticks(x)
ax.set_xticklabels(genes, rotation=90)
ax.legend(loc='upper right')
ax.grid(axis='y', linestyle='--', alpha=0.7)

plt.tight_layout()
plt.savefig("stacked_support_per_gene_mouse.png", dpi=300)
plt.show()
