"""
Version history
- version 1   15-03-2023  K KRISHNA SUDHAKAR
- version 2   09-08-2023 M GRUZIN 
- version 3   18-10-2023 K KRISHNA SUDHAKAR (this version)

GC Profile

    Input:  json target file: SEQ_FILE
            which strand: WHICH_STRAND
            window size: WINDOW

    Output: plot of GC content in window-sized subseq at each position of SEQ
            csv file of GC content in window-sized subseq at each position of SEQ
"""

from curses import window
import os
import json
import numpy
import pandas as pd
import matplotlib.pyplot as plt
from Bio.SeqUtils import gc_fraction
import csv


def open_json_file(json_file, which_strand):
    """
    Open json file and return sequence strand
    Input: json_file, the target file
    Output: sequence strand as string (either forward or reverse strand)

    """
    try:
        #open json file
        data = json.load(open(json_file))

        #initialise sequence variable
        seq = ''

        #check which strand to use
        if which_strand == 'sequence':
            print("Using sequence strand")
            return data['sequence']

        else:
            print("Using sequence_rc strand")
            return data['sequence_rc']

    except Exception as e:
        print(f"Error loading input json file: {e}")
        return None


def gc_profile_calc(input_sequence, window_len, threshold):
    """
    Calculate GC content in window-sized subseq at each position of the sequence
    Input: input_sequence, the sequence strand
    Output: gc_values, a list of GC content in window-sized subseq at each position of the sequence
    """
    try:
        # initialise gc_values list
        gc_values = []

        # initialise high_gc_regions list
        high_gc_regions = []
        region_start = None
        region_end = None

        # loop through sequence to calculate gc content
        for i in range(len(input_sequence) - window_len + 1):
            # extract subsequence (i to i + window_len)
            subseq = input_sequence[i : i + window_len]

            # calculate gc percentage of subsequence
            gc_content_subseq = gc_fraction(subseq)
            gc_perc_subseq = gc_content_subseq * 100

            # append gc content to gc_values list
            gc_values.append(gc_perc_subseq)

            # get start and end of regions that are above threshold
            if region_start is None and gc_perc_subseq >= threshold:
                region_start = i
                region_end = None

            if region_end is None and region_start is not None and gc_perc_subseq < threshold:
                region_end = i   # or end = i + window_len ?
                high_gc_regions.append((region_start, region_end))
                region_start = None

        # print(gc_values)
        # print(f'high_gc_regions: {high_gc_regions}')
        return (gc_values, high_gc_regions)

    except Exception as e:
        print(f"Error calculating GC content: {e}")
        return None


def create_csv_file(seq, list_of_gc_values, window_len, which_strand, offset, circular, file):
    """
    Creates and stores GC values to a CSV file
    Input: list_of_gc_values, a list of GC content in window-sized subseq at each position of the sequence
            window_len, the window size
            which_strand, the strand used
            offset
            circular, true or false
    Output: a csv file with GC values
    """

    csv_filename = file + '.csv'
    csv_file = open(csv_filename, 'w')
    csv_file_writer = csv.writer(csv_file)

    header = ['sub_sequence', 'target_start', 'offset', 'strand', 'window_length', 'GC_%']

    gc_values_csv_list = []
    target_start = 0
    for item in list_of_gc_values:
        if circular and offset > len(seq)-window_len: offset = 0
        gc_value_csv = [
            seq[target_start : target_start + window_len],
            target_start,
            offset,
            which_strand,
            window_len,
            round(item, 1)
        ]
        target_start += 1
        offset += 1
        gc_values_csv_list.append(gc_value_csv)

    csv_file_writer.writerow(header)
    csv_file_writer.writerows(gc_values_csv_list)


def plot_gc_profile(list_of_gc_values, input_sequence, window_len, which_strand, file):
    """
    Plot GC profile of input sequence
    Input: list_of_gc_values, a list of GC content in window-sized subseq at each position of the sequence
           input_sequence, the sequence strand
           window_len, the window size
           which_strand, the strand used
    Output: a plot of GC profile of input sequence

    """

    plt.plot(list_of_gc_values)
    plt.title("GC Profile of '%s' with window size = %i" % (input_sequence[:20] + "...", window_len))
    plt.xlabel(f"{which_strand}")
    plt.ylabel("GC percentage %")
    plt.grid()
    #plt.show()

    plt.savefig(file + ".png", format="png")
    plt.savefig(file + ".pdf", format="pdf")


def cluster_high_gc_regions(high_gc_regions, window_len, threshold, file):
    """
    Get list of start and end of regions that are above threshold. 
    Merge regions that are within 100 bases from each other.
    Add 100 bases to the extreme ends of clusters.
    Input: high_gc_regions, a list of start and ends where GC content is above threshold
            window_len, the window siz
            threshold, the minimum percentage of GC content
    Output: a .... file of ..........
    """
    try:

        filename = file + '.txt'
        cluster_file = open(filename, 'w')
        cluster_file.write("start, end\n")

        clustered_regions = []

        current_region = high_gc_regions[0]

        for region in high_gc_regions[1:]:
            # print(f'current_region[0]: {current_region[0]}')
            # print(f'current_region[1]: {current_region[1]}')
            # print(f'region[0]: {region[0]}')
            # print(f'region[1]: {region[1]}\n')

            if region[0] - current_region[1] <= 100:
                current_region = (current_region[0], region[1])
            elif region[0] - current_region[1] > 100:
                clustered_regions.append((current_region[0], current_region[1]))
                cluster_file.write(f"{current_region[0]}, {current_region[1]}\n")
                current_region = (region[0], region[1])

            # if last element, then append the final cluster
            if region == high_gc_regions[len(high_gc_regions)-1]:
                clustered_regions.append((current_region[0], current_region[1]))
                cluster_file.write(f"{current_region[0]}, {current_region[1]}\n")

        # print(clustered_regions)
        cluster_file.close()
    except Exception as e:
        print(f"Error in clustering high GC regions: {e}")
        return None


def load_seq_file_or_meth(file, which_strand):
    if os.path.isfile(file):
        assert(file.endswith(".json"))
        return pd.DataFrame([{"name": os.path.basename(file), "seq": open_json_file(file, which_strand)}])
    else:
        df = pd.read_csv(file + os.sep + "meta_targets.tsv", sep="\t")
        assert(df.shape[0] > 0)
        assert("seq" in df.columns)
        assert("gene" in df.columns)
        assert("rc_seq" in df.columns)
        col_name = "seq" if which_strand == "sequence" else "rc_seq"
        df = df[["gene", col_name]]
        df.columns = ["name", "seq"]
        return df


def final_gc_profile_run_command(seq_file_or_meth, which_strand, window_len, threshold, circular, offset):
    """
    Final gc profile run command
    Input: all parameters prepared by previous functions
    Output: a plot of GC profile of input sequence
    """
    try:
        df = load_seq_file_or_meth(seq_file_or_meth, which_strand)

        for _, i in df.iterrows():
            name = i["name"]
            seq = i["seq"]

            # confirm sequence length
            print(f'\nLength of target sequence is: {len(seq)}\n')

            # if circular, add first window-sized subseq to end of seq
            if circular is True:
                first_window_subseq = seq[0:window_len]
                # print(f'first_window_subseq: {first_window_subseq}')
                seq += first_window_subseq
                print(f'Circularising the sequence. Length of target sequence is now: {len(seq)}\n')

            # calculate gc content
            print(f'Calculating GC content of {window_len}-bp subsequence at each position of the sequence...\n')
            gc_value_list, high_gc_regions = gc_profile_calc(seq, window_len, threshold)

            # plot gc profile (save in /gc_out folder)
            print(f'Plotting GC profile of {name} at each position of the sequence...')
            out_plot = f"./gc_out/{name}_gc_profile"
            plot_gc_profile(gc_value_list, seq, window_len, which_strand, out_plot)
            print(f'GC profile plot saved as {out_plot}.png and {out_plot}.pdf\n')

            # create csv file of gc profile (save in /gc_out folder)
            print(f'Generating CSV file of GC values of {name} at each position of the sequence...')
            out_csv = f"./gc_out/{name}_gc_values"
            create_csv_file(seq, gc_value_list, window_len, which_strand, offset, circular, out_csv)
            print(f'GC values CSV file saved as {out_csv}.csv')

            # get clustered regions where GC content is above threshold
            print(f'\nClustering high_gc_regions above threshold that are within 100 bases from each other...')
            print(f"threshold: {threshold}")
            out_regions = f"./gc_out/{name}_high_gc_regions"
            cluster_high_gc_regions(high_gc_regions, window_len, threshold, out_regions)
            print(f"Clusters saved to {out_regions}.txt")

    except Exception as e:
        print(f"Error running final gc profile run command: {e}")
        return None