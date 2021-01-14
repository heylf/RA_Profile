import statsmodels.stats.multitest as multi
from scipy import stats
import numpy as np
import time
import sys
import multiprocessing

RNA_VALUE_LIST = dict()
RNA_ID_LIST = dict()
RNA_SEQ_PEAK_COORDS = dict()

# Function to obtain the number of lines in a file.
def get_line_count(filepath):
    file = open(filepath, "r")
    count = 0
    for line in file:
        count += 1
    file.close()
    return count

def reg_potential(atac_peak_coords, rna_seq_peak_coords, max_dist):
    start_a = int(atac_peak_coords[1])
    end_a = int(atac_peak_coords[2])

    start_b = int(rna_seq_peak_coords[1])
    end_b = int(rna_seq_peak_coords[2])

    potential = 0.0

    dist = np.min([np.abs(start_a - start_b), np.abs(end_a - start_b), np.abs(start_a - end_b), np.abs(end_a - end_b)])
    relative_dist = dist / max_dist
    if ( relative_dist <= 1.0 ):
        potential = np.exp(-(0.5 + 4 * relative_dist))

    return potential

def init_peak_analysis(i, atac_data, potential_dict, significance_dict):
    print(i)
    potential = 0.0
    coords_a = atac_data[0]
    chr_a = coords_a[0]
    strand_a = coords_a[3]
    for j in range(0, num_rna_peaks):
        # check chromosome and strandness
        coords_b = RNA_SEQ_PEAK_COORDS[RNA_ID_LIST[j]]
        chr_b = coords_b[0]
        strand_b = coords_b[3]
        corr = 0.0
        pval = 1.0
        if ( strandness == 1 ):
            strand_a = strand_b

        if ( chr_a == chr_b and strand_a == strand_b ):
            corr, pval = stats.spearmanr(atac_data[1], RNA_VALUE_LIST[j])
            if ( np.isnan(corr) ):
                corr = 0.0
            if ( np.isnan(pval) ):
                pval = 1.0
            if ( corr < t ):
                potential += reg_potential(coords_a, coords_b, 100000)
    potential_dict[i] = potential
    significance_dict[i] = (1 - stats.norm.cdf(potential))

# calculate distance
# first check if both of them on same chromosome (yes = continue, no = nan value)
# second check strand (yes = continue, no = nan value)
# min(abs(start_atac - start_rna), abs(end_atac - start_rna)) = distrance

print("[START]")

#atac_peaks_path = "/home/florian/Documents/ATAC_Project/atac_peaks.tabular"
atac_peaks_path = "/home/florian/Documents/ATAC_Project/atac_peaks_subset.tabular"
rna_seq_path = "/home/florian/Documents/ATAC_Project/rna_seq.tabular"

output_path = "/home/florian/Documents/ATAC_Project/"

num_atac_peaks = get_line_count(atac_peaks_path) - 1
num_rna_peaks = get_line_count(rna_seq_path) - 1

atac_peaks_file = open(atac_peaks_path, "r")
headline_x = atac_peaks_file.readline()

rna_seq_file = open(rna_seq_path, "r")
headline_y = rna_seq_file.readline()

#blastn_short_file = open(blastn_short_path, "r")

print("Analysing {} ATAC peaks".format(num_atac_peaks))
print("Analysing {} RNA-Seq peaks".format(num_rna_peaks))

atac_value_list = [ [0.0] ] * num_atac_peaks

atac_id_list = [""] * num_atac_peaks

for i in range(0, num_atac_peaks):
    line_atac = atac_peaks_file.readline()
    x = line_atac.strip("\n").split("\t")
    atac_id_list[i] = x[0]
    atac_value_list[i] = [float(val) for val in x[1:]]

atac_peaks_file.close()

for i in range(0, num_rna_peaks):
    line_rna = rna_seq_file.readline()
    y = line_rna.strip("\n").split("\t")
    RNA_ID_LIST[i] = y[0]
    RNA_VALUE_LIST[i] = [float(val) for val in y[1:]]

rna_seq_file.close()

id_delimiter = "|"

atac_peaks_coords_path = "/home/florian/Documents/ATAC_Project/atac_peaks.bed"
atac_peaks_coords = dict()
atac_coords_to_id = dict()
atac_peaks_coords_file = open(atac_peaks_coords_path, "r")

for line in atac_peaks_coords_file:
    values = line.strip("\n").split("\t")
    if ( values[3] not in atac_peaks_coords ):
        coords = [values[0], values[1], values[2], values[5]]
        atac_peaks_coords[values[3]] = coords
        atac_coords_to_id[id_delimiter.join(coords)] = values[3]

atac_peaks_coords_file.close()

rna_seq_peaks_coords_path = "/home/florian/Documents/ATAC_Project/mm10_annotation.bed"
rna_seq_peaks_coords_file = open(rna_seq_peaks_coords_path, "r")

for line in rna_seq_peaks_coords_file:
    values = line.strip("\n").split("\t")
    if ( values[3] not in atac_peaks_coords ):
        RNA_SEQ_PEAK_COORDS[values[3]] = [values[0], values[1], values[2], values[5]]

rna_seq_peaks_coords_file.close()

# Read blastn file if you want to check also if ATAC-Peaks cover specific transcription factors.
blastn_short_path = "/home/florian/Documents/ATAC_Project/blastn_short.tabular"
blastn_file = open(blastn_short_path, "r")

strandness = 1

num_motif_dict = dict()

for line in blastn_file:
    values = line.strip("\n").split("\t")
    coords = values[0].split(id_delimiter)
    coords = coords[1:]

    if(strandness == 1):
        coords[-1] = "."

    if ( id_delimiter.join(coords) in atac_coords_to_id ):
        id = atac_coords_to_id[id_delimiter.join(coords)]
        if( id not in num_motif_dict ):
            num_motif_dict[id] = 0
        else:
            num_motif_dict[id] += 1

blastn_file.close()

t = 0.05

manager = multiprocessing.Manager()
potential_dict = manager.dict()
significance_dict = manager.dict()

pool = multiprocessing.Pool(4)

start = time.time()
for i in range(0, num_atac_peaks):
    atac_data = [atac_peaks_coords[atac_id_list[i]], atac_value_list[i]]
    try:
        pool.apply_async(init_peak_analysis, args=(i, atac_data, potential_dict, significance_dict))
    except:
        print("[ERROR 1] Failed Analysis")

pool.close()
pool.join()

end = time.time()
print(end - start)

potential_file = open(output_path + "/potential_file.tsv", "w")
potential_file.write("ATAC_Peak\tRegulatory_Potential\tp-val\tcorrected_p-val\tNum_Motifs\n")

corrected_pvals = multi.multipletests(significance_dict.values(), alpha=0.05, method='fdr_bh', is_sorted=False, returnsorted=False)

# List of improtant ATAC_peaks
# all genes with high correlation --> calculate from them the regulatory potential for the atac peak
# pvalue from the cdf of a normal distribution
for i in range(0, len(atac_id_list)):
    peak_id = atac_id_list[i]
    if ( peak_id in num_motif_dict ):
        print(peak_id)
        potential_file.write(("{}\t{}\t{}\t{}\t{}\n").format(peak_id, potential_dict[i], significance_dict[i],
                                                         corrected_pvals[1][i], num_motif_dict[peak_id]))
    else:
        potential_file.write(("{}\t{}\t{}\t{}\t{}\n").format(peak_id, potential_dict[i], significance_dict[i],
                                                         corrected_pvals[1][i], 0))

potential_file.close()
