import pandas as pd
import csv
import sys
from collections import OrderedDict
import pathlib
from pathlib import Path
import re

def find_section_in_file(file, section):
    #Open the relevant csv file
    open_file = csv.reader(open(file))

    # Iterate throw the rows and find the row that reads 'DataType: Units'
    starting_row_found = False
    starting_row = 0
    ending_row = 0
    for index, row in enumerate(open_file):
        print(row[:len(section)])
        if row[:len(section)] == section:
            starting_row = index + 1
            starting_row_found = True

        # Find the first empty row after the section row
        if starting_row_found:
            if row == []:
                ending_row = index - 1
                break
            else:
                if row[0] == '':
                    ending_row = index - 1
                    break

    # Check whether both starting and ending rows were found:
    if starting_row == 0 or ending_row == 0:
        # If not, return None
        return None

    else:
        # If yes, load the median MFI section from the file into a dataframe
        return pd.read_csv(file_path, error_bad_lines=False, skiprows=starting_row, nrows=ending_row-starting_row-1 , header=0)

# Function to find a well in a 96-well plate inside a string
def find_well_in_string(string):
    match = re.search(r'[^A-Z]([A-H][1-9]{1}[0-2]*)[\D]', string)
    if match:
        return match.groups()[0]
    else:
        return False

# Grab a list of folders within the Raw directory, which will be a list of the cell lines
p = Path('../../Data/Raw').glob('*')
cell_lines = [x.name for x in p if not x.is_file()]

# Tack on numbers before the cell lines
cell_lines_str = [f'{i + 1}-{cell_lines[i]}' for i in range(len(cell_lines))]
sep = ", "

# Combine into one string to output to user when asking which cell line they want
cell_lines_output = sep.join(cell_lines_str)

# Ask for which cell line the use wants the data compiled
while True:

    chosen_cell_line = input(f'Choose the number cell line you want the data compiled for:\n '
                             f'{cell_lines_output}:\n')

    # Make sure the use selected a number corresponding to a cell line
    if chosen_cell_line in (str(list(range(1, len(cell_lines) + 1)))):
        cell_line = cell_lines[int(chosen_cell_line) - 1]
        break
    else:
        print("Please choose a number corresponding to one the cell lines\n")

# Grab the csv file for the chosen cell line from the Raw folder
file_path = Path(f'../../Data/Raw/{cell_line}/{cell_line}.csv')

# A dataframe connecting analytes to their bead region ID
analytes_RID_df = find_section_in_file(file=file_path, section=["DataType:", "Units"])

if analytes_RID_df is None:
    sys.exit("Could not find bead region IDs section in file")

analytes_RID_df = analytes_RID_df.drop("Analyte:", axis=1)

# A list of the analytes available
analytes = analytes_RID_df.columns

# This dataframe will associate each sample with its wells
medians_df = find_section_in_file(file = file_path, section=["DataType:", "Median"])

if medians_df is None:
    sys.exit("Could not find medians section in file")

samples_df = medians_df[['Sample', 'Location']]


samples = (samples_df.Sample.unique())

# Grab a list of the raw CSV files associated with the cell line
csv_files_path = Path(f'../../Data/Raw/{cell_line}/CSV/')
csv_files = [x for x in csv_files_path.iterdir() if x.is_file() and x.suffix == '.csv']

# This extracts the wells from the names of the files, allowing for searching for a file associated with a particular well
print(csv_files)
csv_files_wells = [find_well_in_string(str(x)) for x in csv_files]
print(csv_files_wells)
# List where we will put all our dataframes, one dataframe per sample
df_list = []

dataframes = []
for sample in samples:

    # This dictionary will hold the data to be turned into a dataframe
    data_to_save = {'Group': []}

    # Keep track of the number of events captures, important since event numbers are unequal
    event_number = {}

    # Initiate empty series and 0 ev
    for analyte in analytes:
        data_to_save[analyte] = pd.Series([])
        event_number[analyte] = 0
    sample_wells = samples_df['Location'][samples_df['Sample'] == sample].values

    if len(sample_wells) > 0:
        sample_wells = [find_well_in_string(x) for x in sample_wells]
    else:
        sys.exit("Could not match some samples to their repsective wells")

    for well in sample_wells:
        # Grab the data from the file associated with the well
        file = csv_files[(csv_files_wells.index(well))]
        df_file = pd.read_csv(file,skiprows=1)
        for analyte in analytes:
            RID = analytes_RID_df[analyte].values[0]
            events = df_file['RP1'][df_file['RID'] == RID].values
            event_number[analyte] += len(events)

            data_to_save[analyte] = data_to_save[analyte].append(pd.Series(events), ignore_index=True)


    data_to_save['Group'] += [sample.strip()] * max(event_number.values())

    dataframe = pd.DataFrame(data=data_to_save)
    dataframes.append(dataframe)

total_df = pd.concat(dataframes)

total_df.to_csv(Path(f'../../Data/Compiled/{cell_line}.csv'))