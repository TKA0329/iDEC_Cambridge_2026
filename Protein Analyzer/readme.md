# Protein Properties Analyzer
Link to Streamlit Web App: https://proteinanalyzer-tfdxptasphzhddkxctm3ed.streamlit.app/

## How it works

Put your protein sequences in a CSV file with three columns:

1. sequence
2. mutation region (optional, e.g. `3-5, 6, 8, 10`)
3. number of random copies

```csv
sequence,mutation_region,num_random_copies
MPYEKHVEQTVVEKTEQGGSGGSYRHQTEAEAEKIRRELEKQGGGGSGGGGS,3-6,5
RQGGGAPAGGNIGGGQPQGGWGQPQQPQGGNQFSGGAQSRPQ,,1
MKTQRDGHSLGRWSLVLLLLGLVMPLAIIAQVLSYKEAVL,10,8
```

Upload it to the web app and click the mutation/analyze button. The app uses three modes (random, conservative or scan) to create substituted copies per row, then returns a table with properties for each generated sequence: molecular weight, isoelectric point (pI), GRAVY score, aromaticity, instability index, net charge at pH 7, and secondary structure fractions.

## Features

### Mutation modes
Three substitution modes are available when generating sequence variants:

**Random** — each position in the specified region is mutated to a randomly chosen amino acid (any of the 20 standard AAs except the original).

**Conservative** — each position is swapped to a physicochemically similar amino acid only. Polars stay polar, nonpolars stay nonpolar, charged residues stay charged, etc. Based on standard conservative substitution groups.
Note:
In this programme, 
* Cysteine (C) is paired with Serine (S) due to size, but Cysteine is unique because it forms disulphide bridges. If a Cysteine involved in a bridge mutates to Serine, it’s often functionally "non-conservative" because the bridge is lost.

* Proline (P) is paired with Alanine and Glycine. While similar in size, Proline is a "helix-breaker" with a rigid ring. Replacing it often changes the physical backbone of the protein.However, it is included in this programme for informational purposes. 

* Histidine (H) is paired with Arginine (R) and Tyrosine (Y). Histidine is tricky because its charge depends on the local pH, making it an AA that can act as either polar, charged, or aromatic.

**Scan (single-position exhaustive)** — for each position in the region, every possible conservative substitute is generated as a separate variant, changing only that one position at a time. Useful for identifying which specific residue is driving a property change. The `num_copies` column is ignored in this mode.

### Combination counter
An expandable panel shows the total number of possible unique variants for each sequence and region in your CSV, broken down per position. Updates based on whichever mutation mode is selected.


Useful for processing large numbers of sequences at once without having to submit them one by one to tools like ProtParam.

Run locally:

```bash
streamlit run analyzer.py
```
