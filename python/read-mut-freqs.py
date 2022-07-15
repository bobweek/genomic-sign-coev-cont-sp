from curses import meta
import os
import pandas as pd
import numpy as np

# in slim different mutations are tracked each generation
# so the outputted csv with mutations on columns and time on rows
# will have different number of cols for each row
# so we have to deal that uneveness when reading those csv files in

datloc = os.path.expanduser("~/gsccs-data/")

data_file_delimiter = ','

# read in host data
data_file = datloc+"hfs-id.csv"
# The max column count a line in the file could have
largest_column_count = 0
# Loop the data lines
with open(data_file, 'r') as temp_f:
    # Read the lines
    lines = temp_f.readlines()
    for l in lines:
        # Count the column count for the current line
        column_count = len(l.split(data_file_delimiter)) 
        # Set the new most column count
        largest_column_count = column_count if largest_column_count < column_count else largest_column_count
column_names = [i for i in range(0, largest_column_count)]
hfs = np.array(pd.read_csv(datloc+"hfs.csv", delimiter=data_file_delimiter, names=column_names))
hfs_id = np.array(pd.read_csv(datloc+"hfs-id.csv", delimiter=data_file_delimiter, names=column_names))

# read in para data
data_file = datloc+"pfs-id.csv"
# The max column count a line in the file could have
largest_column_count = 0
# Loop the data lines
with open(data_file, 'r') as temp_f:
    # Read the lines
    lines = temp_f.readlines()
    for l in lines:
        # Count the column count for the current line
        column_count = len(l.split(data_file_delimiter)) 
        # Set the new most column count
        largest_column_count = column_count if largest_column_count < column_count else largest_column_count
column_names = [i for i in range(0, largest_column_count)]
pfs = np.array(pd.read_csv(datloc+"pfs.csv", delimiter=data_file_delimiter, names=column_names))
pfs_id = np.array(pd.read_csv(datloc+"pfs-id.csv", delimiter=data_file_delimiter, names=column_names))

# now that we have read the csv files of mutation freqs and ids
# we organize these data into a dataframe that can be read into r
# the dataframe consists of four columns: time, frequency, id, and spp

# in slim each mutation is given a unique id
# so we search for the max id in each spp
# and then iterate through each id starting with one up to the max
# find where that id occurs

tempres = pd.read_csv(datloc+"params.csv")["tempres"][0]

times = []
freqs = []
ids   = []
spp   = []

# organize host freqs
maxid = int(max(hfs_id.flatten()[~np.isnan(hfs_id.flatten())]))
for m in (np.arange(maxid)+1):
    metaloc = np.where(hfs_id==m)
    if np.size(metaloc)>0:
        nfreqs = np.zeros(hfs.shape[0])
        for i in np.arange(np.shape(metaloc)[1]):
            nfreqs[metaloc[0][i]] = hfs[metaloc[0][i],metaloc[1][i]]        
        freqs.append(nfreqs)
        ids.append(np.repeat(m,hfs.shape[0]))
        spp.append(np.repeat(1,hfs.shape[0]))
        times.append(np.arange(hfs.shape[0])*tempres+2)

# organize para freqs
maxid = int(max(pfs_id.flatten()[~np.isnan(pfs_id.flatten())]))
for m in (np.arange(maxid)+1):
    metaloc = np.where(pfs_id==m)
    if np.size(metaloc)>0:
        nfreqs = np.zeros(pfs.shape[0])
        for i in np.arange(np.shape(metaloc)[1]):
            nfreqs[metaloc[0][i]] = pfs[metaloc[0][i],metaloc[1][i]]        
        freqs.append(nfreqs)
        ids.append(np.repeat(m,pfs.shape[0]))
        spp.append(np.repeat(2,pfs.shape[0]))
        times.append(np.arange(pfs.shape[0])*tempres+2)

times = np.array(times).flatten()
spp = np.array(spp).flatten()
ids = np.array(ids).flatten()
freqs = np.array(freqs).flatten()

data = {'t': times,
        'spp': spp,
        'id': ids,
        'frq': freqs}

df = pd.DataFrame(data)

df.to_csv(datloc+'mut-freqs.csv')