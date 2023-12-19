'''
Version 1
08/11/2023 by Kavitha Krishna Sudhakar

This script uses Vienna Suite's RNAlfold to analyse hairpins in the input sequence.
Vienna output is parsed

'''


import sys
import os
import json

from lib.vienna_api_class import ViennaAPI

if __name__ == '__main__':

    # initialise viennaAPI class object
    vienna_api = ViennaAPI()

    SEQ_FILE = sys.argv[1]

    seq_file = open(SEQ_FILE)
    # sequences = json.load(seq_file)
    # seq = sequences['sequence']
    seq = seq_file.readline()

    # print(seq)

    lfold_result = vienna_api.RNALfold(seq, 50)

    vienna_result = lfold_result[0].split('\n')


    del vienna_result[-1]
    del vienna_result[-1]
    del vienna_result[-1]

    # print("======================after deleting=====================")
    # print(vienna_result)

    # create new file to copy lfold results into
    filename = "RNALfold_" + SEQ_FILE + '.lfold'
    lfold_result_file = open(filename, 'w')

    # split each line of lfold result and remove blank spcaes and brackets
    for i in range(len(vienna_result)):
        vienna_split = vienna_result[i].split()
        # print(vienna_split)

        # each vienna_split should have 3 elements
        # as deltaG < '( -9.9)' introduces an unwanted blank space
        if len(vienna_split) == 4:
            if vienna_split[1] == '(':
                del vienna_split[1]

        # print (vienna_split)

        # extract lfold elements of each line 
        lfold_dot_bracket_notation = vienna_split[0]
        lfold_free_energy = vienna_split[1].strip('(').strip(')').strip()
        lfold_start_pos = vienna_split[2]

        # print (lfold_dot_bracket_notation, lfold_free_energy, lfold_start_pos)

        # saving lfold results to file
        lfold_result_file.write('{} {}  {}\n'.format(lfold_dot_bracket_notation, lfold_free_energy, lfold_start_pos))

    
    # sorting lfold results file by highest delta G to lowest delta G
    # os.system('seq_file=SEQ_FILE')
    print('\n=================================== SORTING RNALFOLD LINES =====================================\n')
    # os.system('touch RNALfold_result_sorted.lfold')
    # os.system('cat - RNALfold_STAT3_sequence.txt.lfold | sort -n -k3 > RNALfold_result_sorted.lfold')
    os.system('sort -n -k3 RNALfold_STAT3_plus_GH17J042386.txt.lfold > RNALfold_STAT3_plus_GH17J042386_sorted_by_pos.lfold')
    os.system('head -n 237 RNALfold_STAT3_plus_GH17J042386_sorted_by_pos.lfold | sort -n -k2')
    





