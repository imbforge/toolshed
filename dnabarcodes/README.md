# Calculate pairwise edit distances between DNA barcodes

Small Shiny app to calculate the pairwise Hamming distances for a set of DNA barcodes.
All barcodes need to be of the same length.

Plots a histogram and a heatmap of the pairwise Hamming distances, as well as a sequence logo for the barcode set (using the `ggseqlogo` R package).
The barcode pairs that fall below a desired distance threshold (default = 3) are presented in a table below the plots.

*Input file format*: Tab-separated two-column text file with barcode IDs/names (1st column) and sequences (2nd column).
An example input file is provided in this repository.

