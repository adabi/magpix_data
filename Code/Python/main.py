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
        if row == section:
            starting_row = index + 1
            starting_row_found = True

        # Find the first empty row after the section row
        if starting_row_found:
            if row == []:
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
    match = re.search(r'[A-H]\d{1,2}', string)
    if match:
        return match.group(0)
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
csv_files_wells = [find_well_in_string(str(x)) for x in csv_files]

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





'''


#print(samples_names_df)
data = {"Endothelial": [True, True, True, False, False, False, False, False, False], "Microglial": [False, False, False, True, True, True, False, False, False],
        "Concentration": [0, 1, 20, 0, 1, 20, 0, 1, 20], "Wells":[[] for x in range(9)]}

wells_df = pd.DataFrame(data=data)

for row in samples_names_df.iterrows():
    sample = row[1]["Sample"]
    well = row[1]["Well"]
    split_sample = sample.split(" ")
    if split_sample[1] == "DEP":
        concentration = int(split_sample[2])
        if split_sample[0] == "E":
            (wells_df.Wells[(wells_df["Endothelial"] == True) & (wells_df["Concentration"] == concentration)].values[0]).append(well)
        elif split_sample[0] == "M":
            (wells_df.Wells[(wells_df["Microglial"] == True) & (wells_df["Concentration"] == concentration)].values[0]).append(well)
        else:
            (wells_df.Wells[(wells_df["Microglial"] == False) & (wells_df["Endothelial"] == False) & (wells_df["Concentration"] == concentration)].values[
                 0]).append(well)


print(wells_df)

# Enumerate all files in the 'Data Files' with the csv extension
files = []
for file in os.listdir("Data Files"):
    if file.endswith(".csv"):
        files.append(file)


file_splt = [x.split("_")[2][:-4] for x in files if len(x.split("_"))==3]

full_df = pd.DataFrame(columns=wells_df.columns.tolist()[:-1] + analytes_RID_df.columns.tolist())
print(full_df)
for i, row in wells_df.iterrows():
    well_data = OrderedDict()

    for well in row["Wells"]:
        index = file_splt.index(well)
        file = files[index]
        df = pd.read_csv(os.path.join("Data Files", file), header=1)
        for analyte in analytes_RID_df:
            RID = analytes_RID_df[analyte].values[0]
            medians = df["RP1"][df["RID"] == RID]
            medians.reset_index(drop=True, inplace=True)
            well_data[analyte] = medians
        medians_df = pd.DataFrame(well_data)
        medians_df.insert(0, column="Endothelial", value=[row["Endothelial"] for x in range(len(medians_df))])
        medians_df.insert(1, column="Microglial", value=[row["Microglial"] for x in range(len(medians_df))])
        medians_df.insert(2, column="Concentration", value=[str(row["Concentration"]) + " ug/ml" for x in range(len(medians_df))])

    full_df = pd.concat([full_df, medians_df], ignore_index=True)


full_df.to_csv(path_or_buf="Full Raw Data.csv", index=False)

'''