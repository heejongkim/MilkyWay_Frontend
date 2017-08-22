# MilkyWay - 0.3.9

# Known Issue: 
* In PSM-ID tab, the selected peptide's filtered PSM table for Lorikeet can be woefully presented if the last clicked row's index is bigger than the newly selected peptide's PSM list's total number of rows. Still working on what's the best way to handle this at this point.
* FIX: PSM table's precursor Mz is not M/Z, resulting to mislabel Lorikeet precursor position

# Plan
## FUTURE PLAN
* TO-DO: Integrate Suggestions of M. Choi as in SeeMSstats
* NEW: Visualization of intersecting protein IDs across runs - Interactive UpSet Plot
* NEW: Protein-Protein Interaction Network Visualization

* IMPROVE: Heatmap height auto-adjustment
* IMPROVE: Galaxy Job Submitter Layout improvement (In Progress):
    * DIA or DDA uploader in good condition.
    * DIA+DDA uploader fully functional, but ugly.  Users must be careful to do things in the right order.
    * In the end, the upload tool should be changed away from holding the uploader panels as tabs in a tabBox. Instead, we should populate separate tabs on the dashboard side-panel for uploading... That should allow for proper use of 'conditionalPanel's.
* IMPROVE/FIX: In protein tab, Sequence Coverage Map - 1) Optimize the coverage parsing to show only unique peptides (both un/modificed ones) with PSM counts for each AA as a graph 2) Troubleshoot with multiple lines issue even though there's no overlapping multiple peptides
