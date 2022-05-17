#!/usr/bin/env python

import os
import sys
import pandas as pd
import numpy as np
from tqdm import tqdm
from collections import Counter
import argparse
from datetime import datetime

########################################
## Read in data and set up parameters ##
########################################

#Create command line arguments  
parser = argparse.ArgumentParser(description='Create co-occurrence matrix for EHR code within different time windows.')
parser.add_argument('-i', '--infile', help='Directory for input file (all numeric variables)')
parser.add_argument('-o', '--outfile', help='Directory for output file (code1, code2, count, and window)')
parser.add_argument('-w', '--windows', nargs='+', type=int, help='Enter the time windows (seperate it by whitespace)', default=[0, 1, 6, 13])
parser.add_argument('-c', '--chunk', type=int, help='current chunk number (array number)', default=1)
parser.add_argument('-tc', '--tot_chunks', type=int, help='split data into tot_chuncks of patients to speed up the process', default=1)
args = parser.parse_args()

# window is defined in terms of day difference, e.g. window=0 indicates as the same day
windows = args.windows

# store co-occurence count for each windows
matrices = [Counter() for _ in range(len(windows))]
print(windows)

print('start read in file: ', datetime.today())
# os.getcwd()
# read in data
events = pd.read_csv(args.infile, header=0)

# input chunk and tot_chunks through command line
chunk = args.chunk #chunk = 10
tot_chuncks = args.tot_chunks #tot_chuncks = 4800

# total number of patients and number of patients per chunk
tot_pt = events['PId'].max()
chunk_per_pt = round(tot_pt / tot_chuncks)

# subset data into the chunk from Patient ID 'min_pid' to 'max_pid' 
min_pid = (chunk_per_pt*(chunk-1)+1)
max_pid = chunk*chunk_per_pt

print(min_pid)
print(max_pid)

##################################
## Calculate cooccurence matrix ##
##################################

# subset data to chunks, read the table file and subset based on the calculated ID range
events_perpt = events.loc[(events['PId'] >= min_pid) & (events['PId'] <= max_pid)]

bar = tqdm(total=len(events_perpt))

# calculate cooccurance matrix
for i in range(len(events_perpt)-1): 
    sid, day, code = events_perpt.iloc[i]
    for j in range(i+1, len(events_perpt)):
        nsid, nday, ncode = events_perpt.iloc[j]
        if sid != nsid: break
        diffDay = nday - day   
        if diffDay > windows[-1]: break
        
        #similar to which, looking for which non-overlapping window interval to save the count, e.g. windows=[0,1,6,13], diffDay=3, hence it belongs to the interval [1,6]
        tempLoc=np.where([x>=diffDay for x in windows])[0][0]
        
        #coocurance pair (3,5) and (5,3) should be considered as the same, and we unified them by setting the first code always smaller than the second one
        if code < ncode: matrices[tempLoc][(code,ncode)]+=1
        else: matrices[tempLoc][(ncode, code)]+=1
    bar.update(1)

#########################
## Format output table ##
#########################

# format table to data frame
tempall = [pd.DataFrame.from_dict(matrices[i], orient='index').reset_index() for i in range(len(windows))]
for i in range(len(windows)):
    tempall[i]['window']=windows[i]
    tempall[i]['code1']=[x[0] for x in tempall[i]['index']]
    tempall[i]['code2']=[x[1] for x in tempall[i]['index']]
    
# merge all table into one long table
temp=pd.concat(tempall, axis=0)

# output table dimension
temp.shape

# format column name
temp = temp.rename(columns={0: 'count'})
co_occur = temp.loc[:,['code1','code2','count','window']]


##############################
## Output table to csv file ##
##############################
print('start write out file: ', datetime.today())
outfile = args.outfile
new_outfile = str(outfile.replace('.csv','')+'_'+str(chunk)+'.csv')
co_occur.to_csv(new_outfile, index=False)
