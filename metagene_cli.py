import typer
from typing import List, Optional
import matplotlib
import pandas as pd
from Bio import SeqIO
from matplotlib import pyplot as plt
import numpy as np
import statsmodels.formula.api as smf
import pandas as pd

from scipy import stats
# # from scipy import spatial
# # import json

from collections import Counter
import functools

app = typer.Typer()

matplotlib.rcParams['xtick.labelsize'] = 10
matplotlib.rcParams['ytick.labelsize'] = 10
matplotlib.rcParams['axes.labelsize'] = 10
matplotlib.rcParams['axes.titlesize'] = 10

matplotlib.rcParams['axes.grid'] = True
matplotlib.rcParams['grid.color'] = '0.5'
matplotlib.rcParams['grid.linewidth'] = '0.5'

matplotlib.rcParams['axes.edgecolor'] = '0.25'
matplotlib.rcParams['xtick.color'] = '0'
matplotlib.rcParams['ytick.color'] = '0'

matplotlib.rcParams['xtick.major.width'] = 1
matplotlib.rcParams['ytick.major.width'] = 1
matplotlib.rcParams['ytick.major.size'] = 5
matplotlib.rcParams['xtick.major.size'] = 5
matplotlib.rcParams['axes.spines.right'] = True
matplotlib.rcParams['axes.spines.left'] = True
matplotlib.rcParams['axes.spines.top'] = True
matplotlib.rcParams['axes.spines.bottom'] = True

matplotlib.rcParams['font.family'] = 'DejaVu Sans'
matplotlib.rcParams['font.sans-serif'] = 'Helvetica'
matplotlib.rcParams['font.weight']='normal'
matplotlib.rcParams['axes.axisbelow'] = True


mapping_offset = int(input('what is your mapping offset?'))
utr_length_to_include = int(input('what is your utr length io include?'))
genome_file = input('path to your fasta genome file?')
genome_seq = list(SeqIO.parse(genome_file, 'fasta'))[0]

#filter inputs for later:
filter_length = int(input ('what is your length cutoff (100)?'))
filter_cov = int(input('what is your rpm coverage cutoff (95)?'))
filter_expression = float(input('what is your mean expression cutoff(0.05)?'))


#region_to extend:
region_to_extend = int(input('what is your region to extend (50)?'))


utr_sequence_dict = {}
full_sequence_dict = {}
strand_dict = {}
location_dict = {}
genbank_file = input('path to your genbank file?')
genbank_seq = list(SeqIO.parse(genbank_file, 'genbank'))[0]
for feature in genbank_seq.features:
    if feature.type != 'CDS':
        continue
    elif 'pseudo' in feature.qualifiers:
        continue
    try:
        name = feature.qualifiers['locus_tag'][0] + '_' + feature.qualifiers['gene'][0]
    except KeyError:
        name = feature.qualifiers['locus_tag'][0]
    start = feature.location.start
    end = feature.location.end
    ###Sequence for positive strand genes
    if feature.strand == 1:
        seq = str(genbank_seq.seq[start - utr_length_to_include:end + utr_length_to_include])
    ###Sequence for negative strand genes needs to be reverse complemented
    elif feature.strand == -1:
        seq = str(genbank_seq.seq[start - utr_length_to_include:end + utr_length_to_include].reverse_complement())

    full_sequence_dict[name] = seq

    ###Separately just grab the 5'UTR as its own dictionary because it may be helpful
    # for Api the 3'UTR seems more useful. not sure how to code that - look at next code cell
    utr_seq = seq[:utr_length_to_include]
    utr_sequence_dict[name] = utr_seq

    ###Store strand and locations
    strand_dict[name] = feature.strand
    location_dict[name] = (start, end)

print(len(utr_sequence_dict.keys()))
print(len(full_sequence_dict.keys()))
print(len(strand_dict.keys()))
print(len(location_dict.keys()))




starts = []
stops = []
for name,gene in list(full_sequence_dict.items())[:]:
    starts.append(gene[utr_length_to_include:utr_length_to_include+3])
    stops.append(gene[-3-utr_length_to_include:-utr_length_to_include])
print('####Start codons')
print(Counter(starts))
print('#####Stop codons')
print(Counter(stops))






@app.command()
def start(file_name: Optional[List[str]] = typer.Option(None), file_path_pos: Optional[List[str]] = typer.Option(None), file_path_neg: Optional[List[str]] = typer.Option(None), color: Optional[List[str]] = typer.Option(None)):
    zip_list = zip(file_name, file_path_pos, file_path_neg)
    sample_files = []
    for name, path_pos, path_neg in zip_list:
        sample_files.append((name, path_pos, path_neg))
    sample_names = [i[0] for i in sample_files]

    feature_dict_meta = {}
    for sample_file in sample_files[:]:
        sample_name, fwd, rev = sample_file
        print('##### {}'.format(sample_name))
        feature_dict_meta[sample_name] = {}

        ################################
        ###Essentially get the fwd and rev genome coverage
        fwd_dicty = {}
        rev_dicty = {}
        with open(fwd) as infile:
            for line in enumerate(infile):
                if line[0] > 0:  ###Ignore the first line of the file
                    split_line = line[1].split('\t')
                    fwd_dicty[int(split_line[0]) + mapping_offset] = float(
                        split_line[1])  # Note: mapping offset addition
        print('Done with pos')
        with open(rev) as infile:
            for line in enumerate(infile):
                if line[0] > 0:  ###Ignore the first line of the file
                    split_line = line[1].split('\t')
                    rev_dicty[int(split_line[0]) - mapping_offset] = float(
                        split_line[1])  # Note: mapping offset subtraction
        print('Done with neg')
        for gene_name in full_sequence_dict.keys():
            ####Dealing with positive strand genes first
            if strand_dict[gene_name] == 1:
                ###Get all positions that I care about
                pos = (location_dict[gene_name][0] - utr_length_to_include,
                       location_dict[gene_name][1] + utr_length_to_include)
                if pos[1] < pos[0]:
                    print('found a bug')
                    continue
                sequencing = []
                ###Append read values and if there are none, append zero
                for i in range(pos[0], pos[1]):
                    try:
                        sequencing.append(fwd_dicty[i])
                    except KeyError:
                        sequencing.append(0)
                feature_dict_meta[sample_name][gene_name] = sequencing
            ####And repeat for negative strand genes
            elif strand_dict[gene_name] == -1:
                pos = (location_dict[gene_name][0] - utr_length_to_include,
                       location_dict[gene_name][1] + utr_length_to_include)
                if pos[1] < pos[0]:
                    print('found a bug')
                    continue
                sequencing = []
                for i in range(pos[0] + 1, pos[1] + 1):
                    try:
                        sequencing.append(rev_dicty[i])
                    except KeyError:
                        sequencing.append(0)
                feature_dict_meta[sample_name][gene_name] = sequencing[::-1]  ###Note the CRUCIAL reverse here



    ####Check total read numbers across all of these features (only considering the CDS)
    total_read_dict = {}
    for i in sample_names:
        all_features = []
        for j in feature_dict_meta[i].values():
            all_features.extend(j[utr_length_to_include:-1 * utr_length_to_include])
        print(i, np.sum(all_features))
        total_read_dict[i] = np.sum(all_features)



    ###filters:
    for i in sample_names:
        if i.find('_ribo') == -1:
            continue
        print('##################################')
        print(i, ' inititially has ', len(feature_dict_meta[i]), 'genes (sanity check)')
        to_delete = []
        for j in feature_dict_meta[i]:
            ####Length cutoff is to make sure each CDS is 100 nts long
            if len(feature_dict_meta[i][j]) < filter_length + (utr_length_to_include * 2):
                to_delete.append(j)
            ####Coverage cutoff is to make sure that at least 5% of the values in the CDS are non-zero
            elif np.percentile(feature_dict_meta[i][j][utr_length_to_include:-1 * utr_length_to_include],
                               filter_cov) <= 0:  ###Important parameter here, genes to discard based on coverage
                to_delete.append(j)
            ####Finally cut the lowest expression genes (mean less than 0.5)
            elif np.mean(feature_dict_meta[i][j][utr_length_to_include:-1 * utr_length_to_include]) < filter_expression:
                to_delete.append(j)
        for j in to_delete:
            del feature_dict_meta[i][j]
        print('now has', len(feature_dict_meta[i]), 'after coverage and length threshold')


    tempor_dict = {}
    for sample in sample_names:
        tempor_dict[sample] = feature_dict_meta[sample]

    #tempor_dict = {WT_A_ribo : {gene1;[1,2,3], gene2: [4,5,6]}, TET_A_ribo : {gene1;[1,2,3], gene2: [4,5,6]}}
    genes_in_all_ribo = list(functools.reduce(set.intersection, (set(val) for val in tempor_dict.values())))



    print('##################################')
    print('##################################')
    print('##################################')
    print('There are', len(genes_in_all_ribo), 'genes that appear in all ribosome datasets')

    #####Finally, subset out genes based on operon position to make the meta-analysis more clean

    genes_with_a_preceder = []
    genes_with_a_follower = []
    starts_positive_strand = []
    starts_negative_strand = []
    ends_positive_strand = []
    ends_negative_strand = []
    for gene in strand_dict.keys():
        if strand_dict[gene] == 1:
            starts_positive_strand.append(location_dict[gene][0])
            ends_positive_strand.append(location_dict[gene][1])
            assert starts_positive_strand[-1] < ends_positive_strand[-1]
        if strand_dict[gene] == -1:
            starts_negative_strand.append(location_dict[gene][1])
            ends_negative_strand.append(location_dict[gene][0])
            assert starts_negative_strand[-1] > ends_negative_strand[-1]

    print('Compiled all starts/stops')

    for i, j in enumerate(starts_positive_strand):
        starts_positive_strand[i] = list(range(j - region_to_extend, j + region_to_extend, 1))
    for i, j in enumerate(starts_negative_strand):
        starts_negative_strand[i] = list(range(j - region_to_extend, j + region_to_extend, 1))
    for i, j in enumerate(ends_positive_strand):
        ends_positive_strand[i] = list(range(j - region_to_extend, j + region_to_extend, 1))
    for i, j in enumerate(ends_negative_strand):
        ends_negative_strand[i] = list(range(j - region_to_extend, j + region_to_extend, 1))
    print('Extended all starts/stops by {}'.format(region_to_extend))

    ####Just put everyone together into a bigggg list/s
    all_starts_positive = []
    for i in starts_positive_strand:
        all_starts_positive.extend(i)
    all_starts_negative = []
    for i in starts_negative_strand:
        all_starts_negative.extend(i)

    all_ends_positive = []
    for i in ends_positive_strand:
        all_ends_positive.extend(i)
    all_ends_negative = []
    for i in ends_negative_strand:
        all_ends_negative.extend(i)
    print('Finished concatenating the location of all starts/stops')

    for gene in strand_dict.keys():
        if strand_dict[gene] == 1:
            start, end = location_dict[gene]
            if end in all_starts_positive:
                genes_with_a_follower.append(gene)
            if start in all_ends_positive:
                genes_with_a_preceder.append(gene)
        elif strand_dict[gene] == -1:
            end, start = location_dict[gene]
            if end in all_starts_negative:
                genes_with_a_follower.append(gene)
            if start in all_ends_negative:
                genes_with_a_preceder.append(gene)

    print('Number of genes with a follower', len(genes_with_a_follower))
    print('Number of genes with a preceder', len(genes_with_a_preceder))


    #5' metagene:

    plotting_dict = {}
    full_dict_temp = {}
    winsorize = False
    analysis_type = 'Ribo'
    #####################################
    #####################################
    #####################################
    if analysis_type == 'Ribo':
        samples_to_consider = sample_names[:4]
        genes_to_consider = genes_in_all_ribo
    # elif analysis_type == 'RNA':
    #    samples_to_consider = sample_names[4:]
    #    genes_to_consider = genes_in_all_rna

    for sample in samples_to_consider:
        print(sample)
        temp_array = []
        for gene in genes_to_consider:
            if gene in genes_with_a_preceder:
                continue
            #         if strand_dict[gene] == 1:
            #             continue
            ###
            if winsorize:
                reads = stats.mstats.winsorize(np.array(feature_dict_meta[sample][gene]), axis=0, limits=0.05)
            else:
                reads = feature_dict_meta[sample][gene]
            ###
            meany = np.mean(
                reads[utr_length_to_include:-1 * utr_length_to_include])  ###Calculate the mean based on the CDS only
            temp_array.append([i / meany for i in reads[:utr_length_to_include * 2]])
        print(len(temp_array))
        y_vals = np.mean(np.array(temp_array), axis=0)
        plotting_dict[sample] = y_vals
        full_dict_temp[sample] = temp_array

    print(len(y_vals))

    x_vals = np.arange(-1 * utr_length_to_include, utr_length_to_include, 1)
    fig, ax = plt.subplots(figsize=(8, 6))
    counter = 0
    for i in range(0, len(sample_names)):
        ax.plot(x_vals, plotting_dict[samples_to_consider[counter]], c=color[counter], alpha=0.8, label=sample_names[counter])
        counter += 1
    ax.legend(fontsize=16)
    ax.set_title('{} genes appearing in all 4 {} datasets'.format(len(temp_array), analysis_type), fontsize=20)
    ax.tick_params(labelsize=16)
    ax.set_xlabel('Position relative to start', fontsize=20)
    ax.set_ylabel('Avg normalized reads', fontsize=20)
    maxy = max([max(plotting_dict[i]) for i in samples_to_consider])

    # ax.plot([0, 0], [0, maxy + 0.2], 'k-')
    ax.set_ylim(0, maxy + 0.2)
    ax.grid(False)
    plt.tight_layout()
    if winsorize:
        plt.savefig('winsor.pdf', bbox_inches='tight')
    else:
        plt.savefig('5prime_meta.pdf', bbox_inches='tight')

if __name__ == '__main__':
    app()