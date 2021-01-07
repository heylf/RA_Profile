import statsmodels.stats.multitest as multi
from scipy import stats
import numpy as np
import time

# Function to obtain the number of lines in a file.
def get_line_count(filepath):
    file = open(filepath, "r")
    count = 0
    for line in file:
        count += 1
    file.close()
    return count

def reg_potential(pval, t, atac_peak_coords, rna_seq_peak_coords, strandness, max_dist):
    start_a = int(atac_peak_coords[1])
    end_a = int(atac_peak_coords[2])

    start_b = int(rna_seq_peak_coords[1])
    end_b = int(rna_seq_peak_coords[2])

    potential = 0.0

    if (corr < t ):
        dist = np.min([np.abs(start_a - start_b), np.abs(end_a - start_b)])
        relative_dist = dist / max_dist
        if ( relative_dist <= 1.0 ):
            potential = np.exp(-(0.5 + 4 * relative_dist))

    return potential

# calculate distance
# first check if both of them on same chromosome (yes = continue, no = nan value)
# second check strand (yes = continue, no = nan value)
# min(abs(start_atac - start_rna), abs(end_atac - start_rna)) = distrance

print("[START]")

atac_peaks_path = "/home/florian/Documents/ATAC_Project/atac_peaks.tabular"
rna_seq_path = "/home/florian/Documents/ATAC_Project/rna_seq.tabular"
#blastn_short_path = "/home/florian/Documents/ATAC_Project/blastn_short.tabular"

id_delimiter = "|"

output_path = "/home/florian/Documents/ATAC_Project/"

num_atac_peaks = get_line_count(atac_peaks_path) - 1
num_rna_peaks = get_line_count(rna_seq_path) - 1

atac_peaks_file = open(atac_peaks_path, "r")
headline_x = atac_peaks_file.readline()

rna_seq_file = open(rna_seq_path, "r")
headline_y = rna_seq_file.readline()

#blastn_short_file = open(blastn_short_path, "r")

print(num_atac_peaks)
print(num_rna_peaks)

atac_value_list = [ [0.0] ] * num_atac_peaks
rna_value_list = [ [0.0] ] * num_rna_peaks

atac_id_list = [""] * num_atac_peaks
rna_id_list = [""] * num_rna_peaks

for i in range(0, num_atac_peaks):
    line_atac = atac_peaks_file.readline()
    x = line_atac.strip("\n").split("\t")
    atac_id_list[i] = x[0]
    atac_value_list[i] = [float(val) for val in x[1:]]

atac_peaks_file.close()

for i in range(0, num_rna_peaks):
    line_rna = rna_seq_file.readline()
    y = line_rna.strip("\n").split("\t")
    rna_id_list[i] = y[0]
    rna_value_list[i] = [float(val) for val in y[1:]]

rna_seq_file.close()

atac_peaks_coords_path = "/home/florian/Documents/ATAC_Project/atac_peaks.bed"
atac_peaks_coords = dict()
atac_peaks_coords_file = open(atac_peaks_coords_path, "r")

for line in atac_peaks_coords_file:
    values = line.strip("\n").split("\t")
    if ( values[3] not in atac_peaks_coords ):
        atac_peaks_coords[values[3]] = [values[0], values[1], values[2], values[5]]

atac_peaks_coords_file.close()

rna_seq_peaks_coords_path = "/home/florian/Documents/ATAC_Project/mm10_annotation.bed"
rna_seq_peaks_coords = dict()
rna_seq_peaks_coords_file = open(rna_seq_peaks_coords_path, "r")

for line in rna_seq_peaks_coords_file:
    values = line.strip("\n").split("\t")
    if ( values[3] not in atac_peaks_coords ):
        rna_seq_peaks_coords[values[3]] = [values[0], values[1], values[2], values[5]]

rna_seq_peaks_coords_file.close()

correlation_file = open(output_path + "/correlation_file.tsv", "w")
correlation_file.write("ATAC_Peak\t")
correlation_file.write("\t".join(rna_id_list) + "\t")

pval_file = open(output_path + "/pval_file.tsv", "w")
pval_file.write("ATAC_Peak\t")
pval_file.write("\t".join(rna_id_list) + "\t")

strandness = 1

potential_list = [-1] * len(atac_id_list)
significance_list = [-1] * len(atac_id_list)

for i in range(0, num_atac_peaks):
    potential = 0.0
    correlation_file.write(atac_id_list[i] + "\t")
    pval_file.write(atac_id_list[i] + "\t")
    coords_a = atac_peaks_coords[atac_id_list[i]]
    chr_a = coords_a[0]
    strand_a = coords_a[3]
    start = time.time()
    for j in range(0, num_rna_peaks):
        # check chromosome and strandness
        coords_b = rna_seq_peaks_coords[rna_id_list[j]]
        chr_b = coords_b[0]
        strand_b = coords_b[3]
        corr = 0.0
        pval = 1.0
        if ( strandness == 1 ):
            strand_a = strand_b

        if ( chr_a == chr_b and strand_a == strand_b ):
            corr, pval = stats.spearmanr(atac_value_list[i], rna_value_list[j])
            if ( np.isnan(corr) ):
                corr = 0.0
            if (np.isnan(pval)):
                pval = 1.0
            potential += reg_potential(pval, 0.05, atac_peaks_coords[atac_id_list[i]],
                                       rna_seq_peaks_coords[rna_id_list[j]], 1, 100000)
        correlation_file.write(str(corr) + "\t")
        pval_file.write(str(pval) + "\t")
    significance = (1-stats.norm.cdf(potential))/stats.norm.cdf(0)
    potential_list[i] = potential
    significance_list[i] = significance
    end = time.time()
    print(end - start)


potential_file = open(output_path + "/potential_file.tsv", "w")
potential_file.write("ATAC_Peak\tRegulatory_Potential\tp-val\n")

corrects_pvals = multi.multipletests(significance_list, alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

for id in range(0, len(atac_id_list)):
    potential_file.write(("{}\t{}\t{}\t{}\n").format(atac_id_list[i], potential_list[i], significance_list[i], corrects_pvals[i]))

potential_file.close()
pval_file.close()
correlation_file.close()

# List of ATAC_peaks with corr values to gene
# List of ATAC_peaks with pvalue to gene

# List of improtant ATAC_peaks
# all genes with high correlation --> calculate from them the regulatory potential for the atac peak
# pvalue from Walds-test

# List of important genes
# all atac_peaks with high corllation --> calculate from them the regulatory potential for the gene
# pvalue from Walds-test