"""
This is the top 'wrapper' script for the heatmap module
All user interaction occurs here, to protect the code from unintended meddling
Usage: 'python3 run_heatmap.py <experiment name>'
    # nb: if using Docker, instead use: 
    docker run -v ${PWD}:/src --rm -t heatmap python3 run_heatmap.py <data>

"""

import os
import sys
import time

sys.path.append('./src/') # This is so we have access to the other scripts
from heatmap import final_heatmap_run_command

############################################################
# Adjust input parameters here #############################
############################################################

# probes and target
PROBES_NAME = "SHH_transcript_variant_1"
TARGET_NAME = "ZRS_900bp"

# heatmap labels
Y_AXIS_LAB = "ZRS-900bp-cis(minus)"
Y_AXIS_TICKS = 20

# heatmap col name
# Number of items in COL_NAMES must agree with number of csv files in INPUT_DIR_NAME
# nb: make sure this is one chr string ("a,b,c"), instead of multiple str's e.g. ("a","b","c")
# COL_NAMES = "1, 2, 3" # SV40_short_test_file2"
# todo automate

# dummy/test data files (only relevant if running code from IDE not from terminal)
DEBUG_MODE = 'Terminal'
# DEBUG_MODE = 'Interactive' # use this to run from IDE like VSCode
# For 'Interactive' mode, test data is stored in this directory
INPUT_DIR_NAME = './data'

###########################################################
# CODE STARTS HERE ########################################
###########################################################

if __name__ == "__main__":
    START_TIME_TOTAL = time.time()

    MODULE_NAME = "heatmap.py"

    # Where to find the input CSV files?
    len_args = len(sys.argv)
    for i in range(len_args):
        print(f'i: {i}; sys.argv[{i}]: {sys.argv[i]}')

    # Normal mode is to run via Terminal
    if DEBUG_MODE == 'Terminal':
        INPUT_DIR_NAME = sys.argv[1]
    elif DEBUG_MODE == 'Interactive':
        pass
    print(f'Debug: running in {DEBUG_MODE} mode; will use test data in {INPUT_DIR_NAME} directory')

    if not os.path.exists(INPUT_DIR_NAME):
        raise FileNotFoundError("Can't find INPUT_DIR_NAME (argv[1]): {INPUT_DIR_NAME}")

    os.system("mkdir -p out")

    (csv_file_name, r_script_name, heatmap_pdf_name) = \
        final_heatmap_run_command(
        MODULE_NAME,
        PROBES_NAME,
        TARGET_NAME,
        INPUT_DIR_NAME,
        Y_AXIS_LAB,
        Y_AXIS_TICKS
    )

    check_pdf_size = os.path.getsize(heatmap_pdf_name)
    if check_pdf_size < 10000 :
        print(f'Expecting larger PDF size (PDF is {check_pdf_size/1000} KB). Check RScript in "/out" for any issues.')
    else :
        print('Successfuly generated heatmaps and stored in "/out".')

    END_TIME_TOTAL = time.time()
    ELAPSED_TIME_TOTAL = END_TIME_TOTAL - START_TIME_TOTAL
    print(f'\n>>> Total Elapsed time for {MODULE_NAME}: {ELAPSED_TIME_TOTAL:.4f} sec ({ELAPSED_TIME_TOTAL:.2E})')
    print(40 * '-' + '\n')
