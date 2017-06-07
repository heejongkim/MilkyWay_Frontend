# MilkyWay - 0.3.3

# Known Issue: 
* In PSM-ID tab, the selected peptide's filtered PSM table for Lorikeet can be woefully presented if the last clicked row's index is bigger than the newly selected peptide's PSM list's total number of rows. Still working on what's the best way to handle this at this point.
* FIX: PSM table's precursor Mz is not M/Z, resulting to mislabel Lorikeet precursor position

# Plan
## FUTURE PLAN

* NEW: Visualization of intersecting protein IDs across runs - Interactive UpSet Plot
* NEW: Protein-Protein Interaction Network Visualization

* IMPROVE: Heatmap height auto-adjustment
* IMPROVE: Galaxy Job Submitter Layout improvement
* IMPROVE/FIX: In protein tab, Sequence Coverage Map - 1) Optimize the coverage parsing to show only unique peptides (both un/modificed ones) with PSM counts for each AA as a graph 2) Troubleshoot with multiple lines issue even though there's no overlapping multiple peptides
