import re
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Path to the input file containing CIGAR strings
input_file = "/Users/milicaurosevic/Desktop/TP53_test/cigar_SRR8551562.txt"
output_csv = "/Users/milicaurosevic/Desktop/TP53_test/max_insertions_filtered50.csv"
output_plot = "/Users/milicaurosevic/Desktop/TP53_test/histogram_insertions_filtered50.png"

# List to store results
results = []

# Open the file and process each CIGAR string
with open(input_file, "r") as f:
    for idx, line in enumerate(f, start=1):  # Index starts from 1 to number each read
        cigar = line.strip()

        # Find all insertions (I) in the CIGAR string
        insertions = [int(x) for x in re.findall(r"(\d+)I", cigar)]
        
        # Find the maximum insertion, or 0 if none exist
        max_insertion = max(insertions, default=0)

        # Keep only insertions less than 50 nt
        if max_insertion < 50:
            results.append((idx, max_insertion))

# Create a DataFrame to store results
df = pd.DataFrame(results, columns=["Read Number", "Max Insertion"])

# Save the filtered table to a CSV file
df.to_csv(output_csv, index=False)
print(f"Filtered results saved to: {output_csv}")

# Define bin edges with increments of 1 (from 0 to 50)
bins = np.arange(0, 51, 1)  # Step of 1, up to slightly beyond 50

# Create a clean histogram
plt.figure(figsize=(10, 6))
plt.hist(df["Max Insertion"], bins=bins, alpha=0.7)  # Removed edgecolor="black"

# Set x-axis ticks at intervals of 5
plt.xticks(np.arange(0, 51, 5))

# Add labels
plt.xlabel("Longest insertion per read (< 50 nt)")
plt.ylabel("Number of reads")
plt.title("Distribution of longest insertions per read (filtered)")

# Remove grid lines for a cleaner look
plt.grid(False)

# Save histogram to a file
plt.savefig(output_plot, dpi=300)
print(f"Filtered histogram saved to: {output_plot}")

# Show histogram
plt.show()
