import pandas as pd
import os
import csv
import sys
from collections import OrderedDict
from pathlib import Path
#open_file = csv.reader(open(os.path.join("Data Files", "Rat Cyto 2.csv")))
p = Path('../../Data/Raw').glob('*')
cell_lines = [x for x in p if not x.is_file()]

print(cell_lines)

'''
# Iterate throw the rows and find the row that reads 'DataType: Units'
starting_row_found = False
starting_row = 0
ending_row = 0
for index, row in enumerate(open_file):
    if row == ['DataType:', 'Units']:
        starting_row = index + 1
        starting_row_found = True
    # Find the first empty row after the row that reads 'DataType: Median'
    if starting_row_found:
        if row == []:
            ending_row = index - 1
            break

# Check whether both starting and ending rows were found:
if starting_row == 0 or ending_row == 0:
    # If not, throw an error and exit
    sys.exit("No medians section was found in the compiled data file.")

else:
    # If yes, load the median MFI section from the file into a dataframe
    analytes_RID_df = pd.read_csv(os.path.join("Data Files", "Rat Cyto 2.csv"), error_bad_lines=False, skiprows=starting_row, nrows=ending_row-starting_row-1 , header=0)

analytes_RID_df = analytes_RID_df.drop("Analyte:", axis=1)
print(analytes_RID_df)


#print(analytes_RID_df[analytes_RID_df.loc["Analyte:"] == "G-CSF"])
samples_names_df = pd.read_csv(os.path.join("Data Files", "samples.csv"))
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