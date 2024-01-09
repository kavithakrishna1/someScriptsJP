"""
Pre-processing src code for generating heatmaps

Performs pre-processing in 3 functions:
FUNCTION 1: find_target_file (return: len_target)
1. Determines length of target (
    input: directory/file path and module name
    module: opens json target file and determines length of target
    output: list of target file and its length (list_of_file_and_length)
    )

FUNCTION 2: run_prepare_heatmap_files (return: probe_array)
2. Determines length of probes (
    input: directory/file path and module name
    module: opens csv files in directory and determines length of file, then subtracts this value from len target (from 1.)
    output: tuple of probe columns and their length (list_of_file_and_length)
    )
3. Creates new csv file to extract data for heatmaps
    input: directory/file path and module name
    module:
        Sorts the files in the input directory
        Opens each csv file
        Extracts the relevant ΔG column
        Writes that ΔG column into a new csv file
    output: new CSV file containing the relevant ΔG columns from each separate BoAB Finder file

FUNCTION 3: final_heatmap_run_command
4. Final function, which is called from run_heatmap.py code
    input: returns from previous 3 functions and other inputs specified in run_heatmaps.py
    output: R, csv and heatmap pdf files


"""

import os
import csv
import json
import pathlib
import re


def sorted_nicely(alpha_num_list):
    """
    Implementes a "natural sort":
    filenames sort the way a person would expect them to sort:
    - 1, 11, 100 will sort in that order
    - and not as 1, 100, 11
    taken from stackoverflow.com (https://stackoverflow.com/questions/2669059/how-to-sort-alpha-numeric-set-in-python)
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key)]
    return sorted(alpha_num_list, key = alphanum_key)


def determine_target_length_and_seq(
        file
    ):
    """
    Opens target file in json format and determines target length and sequence

    Input: file is your chosen target file in json format
    Returns: list containing tuples (len_target, target sequence)
    """
    with open(file, 'r') as local_file:
        json_data = json.load(local_file)
    len_target = len(json_data["sequence"])
    target_sequence = json_data["sequence"]
    print(f'Target sequence is: {target_sequence}')
    return (len_target, target_sequence)


def find_target_file(
    dir_name
    ):
    """
    Iterates through a directory, and determines length of target(s) from the json file(s)

    Input: dir_name is directory holding target files in json format
    Returns: list containing tuples(file_name, len_target)
    """
    print(f'\nStep 1: Determining the length of target used from json target file in {dir_name}')
    assert os.path.exists(dir_name), f"Can't find directory of input files in {dir_name}"

    # initialise list_of_file_and_length to hold the info to be returned
    list_of_file_and_length = []

    # iterative search through directory, looking for .json files with search string 'target'
    # for which_file in sorted(os.listdir(dir_name)):
    for which_file in sorted_nicely(os.listdir(dir_name)):
        # print(f'Debug: processing {which_file}')

        # skip hidden files (startswith ".")
        if not which_file.startswith('.'):

            # skip non-.json files
            file_extension = pathlib.Path(which_file).suffix
            if file_extension != '.json':
                continue

            search_string = "target"
            if search_string not in which_file:
                continue

            # joining path to easily identify it
            full_json_file_name = os.path.join(dir_name, which_file)
            assert os.path.exists(full_json_file_name), print(f'Cannot find file: {full_json_file_name}')

            len_target = determine_target_length_and_seq(full_json_file_name)[0]
            pair = [full_json_file_name, len_target]
            list_of_file_and_length.append(pair)

            # temp print statement (could be removed from production version?)
            print(f"The length of the target ({which_file}) is: {len_target}.")
            print("Storing file name and length as variable, 'list_of_file_and_length'")

    # finished - make sure there was some work actually done
    if not list_of_file_and_length:
        raise FileNotFoundError (f'Cound not find any eligible json files in directory {dir_name}')

    if len(list_of_file_and_length) > 1:
        raise ValueError (
            f'There are MULTIPLE eligible json files in directory {dir_name}\
            Please choose 1 target file only'
        )

    return len_target


def transpose_array(
        data_array
    ):
    """
    Function to transpose columns->rows, rows->columns

    Input: data_array as columns
    Returns: transposed_array which has been transposed to rows
    """

    print(f'Reached transpose_array; Data_array is size {len(data_array)}')

    # (preparatory step): set max_row and max_col to 0 to ensure dim > 0,0
    max_row = 0
    max_col = 0

    # determine the dimensions of input matrix
    for nrow, row in enumerate(data_array):
        if nrow > max_row:
            max_row = nrow

        for ncol, col in enumerate(row):
            if ncol > max_col:
                max_col = ncol

    print(f'max_row: {max_row}; max_col: {max_col}')

    # build output matrix
    transposed_array = []

    for ncol in range(max_col + 1):
        transposed_array.append([])
        for nrow in range(max_row + 1):
            transposed_array[ncol].append("NA")

    # now populate the matrix
    for nrow, row in enumerate(data_array):
        for ncol, col in enumerate(row):
            transposed_array[ncol][nrow] = data_array[nrow][ncol]

    return transposed_array


def save_transposed_data_array_to_csv(transposed_array, out):
    """
    Receives transposed_array
    Opens a csv file in the default output directory ./out/
    Writes the data array into that file
    """

    with open(out, 'w') as w:
        csv_file_writer = csv.writer(w)
        csv_file_writer.writerows(transposed_array)


def prepare_heatmap(input_dir_name, len_target, out, target_from, target_to):
    """
    Iterates through a directory, finds csv files from a boab_finder analysis
    Opens boab_finder json files, extracts the ΔG of all alignments
    Returns these ΔG, plus also the length of the (calculated) probe

    Input: 3 inputs
    - input_dir_name: location of the json files
    - len_target: the (fixed) length of the target region
    - out: the name of the csv file holding the detla G values for each experiment

    Returns: list containing tuple (data_array, probe_array)
    - data_list: list of ΔGs from multiple probe vs target boab_finder run
    - probe_array: list of length of probes used
    """

    # initialise data array. This will hold the lists of delta_G values from the csv files
    data_array = []  # This will hold the ΔG values
    probe_array = []  # This will hold the (file_name, probe_length) list of tuples

    min_target_start = float("inf")
    max_target_end = 0

    # access input files in input directory
    print(f'\nStep 2: Determining the length of the probe in each csv file in {input_dir_name}')
    assert os.path.exists(input_dir_name), f"Can't find directory of input files {input_dir_name}"

    # For each finder algorithm's output in the input directory
    for which_file in sorted_nicely(os.listdir(input_dir_name)):
        # Skip hidden files (startswith ".")
        if not which_file.startswith('.'):
            # Skip non-csv files
            file_extension = pathlib.Path(which_file).suffix
            if file_extension != '.csv':
                continue

            # Remaining files must be csv files
            full_input_csv_file = os.path.join(input_dir_name, which_file)
            assert os.path.exists(full_input_csv_file), print(f'Cannot find file: {full_input_csv_file}')

            # Initialise data_list to store the delta G scores for each csv file
            file_name = pathlib.Path(which_file).stem
            data_list = [file_name]

            # Open the csv file and read it
            with open(full_input_csv_file, newline='') as csv_file:
                csv_file_contents = csv.reader(csv_file, quotechar='"')

                # Now go through the file, row by row, and extract the ΔG from column 13
                for i_row, row in enumerate(csv_file_contents):
                    nrows_in_file = i_row
                    
                    # Extract ΔG for duplex structure. This is located in column 14 (Python n=13)
                    delta_G_duplex = row[13]

                    if i_row == 0:
                        assert delta_G_duplex == 'duplex deltaG', print(
                            f'Missing header in row" {nrows_in_file}"'
                            )
                    else:
                        data_list.append(delta_G_duplex)

                        target_start = int(row[8])
                        target_end = int(row[9])

                        # Outside the boundary?
                        if target_start < target_from or target_end > target_to:
                            continue

                        min_target_start = min(min_target_start, target_start)
                        max_target_end = max(max_target_end, target_end)

                # Calculate the length of the probe
                len_probe = abs(len_target - nrows_in_file)

                # Add to probe_array
                new_entry_for_list = (file_name, len_probe)
                probe_array.append(new_entry_for_list)

            data_array.append(data_list)
        else:
            # skip hidden files
            print(f'found a dot_file: {which_file}; skipping...')

    size_of_array = len(data_array)
    print(f'array size is {size_of_array}')

    print('\nFinal results returned in probe_array are:')
    for nitem, item in enumerate(probe_array):
        print(f'nitem: {nitem+1} : file_name: {item[0]}, probe length: {item[1]}')

    # Transpose data array
    print('\nStep 3: Generating an input file for heatmaps from several BoAB Finder output csv files')
    assert len(data_list) > 0, print('There is no data to transpose or process')
    transposed_array = transpose_array(data_array)

    # Write data array into output file
    save_transposed_data_array_to_csv(transposed_array, out)

    return (data_array, probe_array, min_target_start, max_target_end)


def probe_vector(input_probe_variable):
    """
    Builds the probe_vector
    """
    probe_vector_variable = [item[1] for item in input_probe_variable]
    print(f'probe vector is: {probe_vector_variable}')
    return probe_vector_variable


def final_heatmap_run_command(
    module_name,
    probes_name,
    target_name,
    input_directory_name,
    y_axis,
    y_axis_ticks,
    target_from = 0,
    target_to = float("inf")):
    """
    The actual module where the R script is edited, constructed and executed
    """
    print(f'=== {module_name} ====')
    assert os.path.exists(input_directory_name), print(f'Unable to find input directory {input_directory_name}')

    # Setting up experiment
    experiment_name = f'{probes_name}_{target_name}'
    title_heatmap_each_probe_mean = f"Heatmap visualising BoABs between {probes_name} \n \
        and {target_name} (mean of each probe)"
    title_heatmap_whole_matrix_mean = f"Heatmap visualising BoABs between {probes_name} \n \
        and {target_name} (mean of entire dataset)"

    os.system("mkdir -p out")
    TMP_CSV_FILE = f'out/heatmap_csv_input_{experiment_name}.csv'

    TARGET_LEN = find_target_file(input_directory_name)

    # Generate a CSV file for running heatmap in R
    (data_array, probe_array, min_target_start, max_target_end) = prepare_heatmap(
        input_directory_name,
        TARGET_LEN,
        TMP_CSV_FILE,
        target_from,
        target_to)

    PROBE_VECTOR_AS_STR = str(probe_vector(probe_array))[1:-1]
    print("\nStep 4: Creating Heatmaps")
    assert os.path.exists(TMP_CSV_FILE), print(f'Cannot find TMP_CSV_FILE {TMP_CSV_FILE}')

    # Where the temporary R-scipt file will goto
    TMP_R_FILE = f'out/heatmap_script_{experiment_name}.R'

    # Generate the heatmap visualization in R
    pdf_name = f'out/heatmap_pdf_output_{experiment_name}.pdf'

    col_names = []
    for i in range(len(data_array)):
        col_names.append(str(i + 1))
    col_names = ", ".join(col_names)
    
    # Estimated y-axis increments to the specified number of ticks
    y_axis_increment = int((max_target_end - min_target_start) / y_axis_ticks)

    with open("./src/rscript_automated_heatmap.R") as r:
        code = r.read()
        assert "<<INPUT_CSV>>" in code, print(
            'Cannot find <<INPUT_CSV>> placeholder in R template <./src/230629mg_heatmap_automation_v6.R>'
            )
        code = code.replace("<<INPUT_CSV>>", TMP_CSV_FILE)
        #assert os.path.exists(TMP_CSV_FILE), print(f'Cannot find TMP_CSV_FILE {TMP_CSV_FILE}')
        code = code.replace("<<INPUT_LEN_TARGET>>", str(TARGET_LEN))
        code = code.replace("<<PROBE_VECTOR_FOR_R>>", PROBE_VECTOR_AS_STR)
        code = code.replace("<<PROBE_NAME>>", probes_name)
        code = code.replace("<<NAME_PDF_OUTPUT>>", pdf_name)
        code = code.replace("<<INPUT_TARGET_START>>", str(min_target_start))
        code = code.replace("<<INPUT_TARGET_END>>", str(max_target_end))
        code = code.replace("<<INPUT_TITLE_HEATMAP_EACH_PROBE_MEAN>>", title_heatmap_each_probe_mean)
        code = code.replace("<<INPUT_TITLE_HEATMAP_MATRIX_MEAN>>", title_heatmap_whole_matrix_mean)
        code = code.replace("<<INPUT_Y_AXIS_LAB>>", y_axis)
        code = code.replace("<<INPUT_Y_AXIS_INCREMENT>>", str(y_axis_increment))
        code = code.replace("<<COLNAMES>>", col_names)
        with open(TMP_R_FILE, "w") as w:
            w.write(code)

    # Run the R script
    os.system(f"Rscript {TMP_R_FILE}")

    return (TMP_CSV_FILE, TMP_R_FILE, pdf_name)
