import matplotlib.pyplot as plt

# Raw OD600 readings (replicates), not blank-subtracted
raw = {
    "pBAD":   [[0.062, 0.067, 0.064], [0.079, 0.079], [0.142, 0.154, 0.154], [0.155, 0.159, 0.158]],
    "Ori":    [[0.068, 0.064, 0.067], [0.079, 0.081], [0.147, 0.141, 0.151], [0.162, 0.172, 0.163]],
    "ssDNA":  [[0.061, 0.065, 0.071], [0.068, 0.069], [0.073, 0.075, 0.077], [0.085, 0.086, 0.085]],
}

blank = 0.041
conversion_factor = 4.1  # plate reader -> spectrophotometer-equivalent

time_labels = ["T0", "T1", "T2", "T3"]
time_minutes = [0, 46, 116, 181]  # elapsed time from T0, for reference/annotation

# Compute blank-subtracted, spectrophotometer-equivalent means per time point
processed = {}
for construct, timepoints in raw.items():
    means = []
    for reps in timepoints:
        mean_raw = sum(reps) / len(reps)
        mean_corrected = (mean_raw - blank) * conversion_factor
        means.append(mean_corrected)
    processed[construct] = means

# Plot
fig, ax = plt.subplots(figsize=(6.5, 4.5))

colors = {"pBAD": "#B3541E", "Ori": "#2A6F97", "ssDNA": "#5E8C3A"}
labels = {"pBAD": "pBAD", "Ori": "WT", "ssDNA": "Library (ssDNA)"}
markers = {"pBAD": "s", "Ori": "o", "ssDNA": "^"}

for construct, values in processed.items():
    ax.plot(time_labels, values, marker=markers[construct],
             label=labels[construct], color=colors[construct])

ax.set_xlabel("Time point")
ax.set_ylabel(r"OD$_{600}$")
ax.set_title("Growth of pBAD, WT and library cultures (T0-T3)")
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig("od_growth_curve.png", dpi=300, bbox_inches="tight")
plt.show()

# Print processed values for reference / reporting
print(f"{'Timepoint':<10}{'pBAD':<10}{'Ori':<10}{'ssDNA':<10}")
for i, t in enumerate(time_labels):
    print(f"{t:<10}{processed['pBAD'][i]:<10.3f}{processed['Ori'][i]:<10.3f}{processed['ssDNA'][i]:<10.3f}")