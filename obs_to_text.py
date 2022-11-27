#!/usr/bin/python                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
import os
import sys
import numpy as np
import argparse

parser = argparse.ArgumentParser(description='Pulsar profile alignment for variability studies')
parser.add_argument('-p','--pulsar', help='Name as pulsar', required=True)
parser.add_argument('-d','--directory', help='Directory in which the observation files are located', required=True)
parser.add_argument('-b','--phasebins', help='Number of phase bins in the observations', type=int, required='True')
parser.add_argument('-f','--obsfreq', help='Approximate observation frequency (MHz)', type=int, required='True')

# PASSING THE ARGUMENTS TO VARIABLES                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               

args = parser.parse_args()
pulsar = args.pulsar
directory = args.directory
bins = args.phasebins
obs_freq = args.obsfreq

# ALLOW GREP IN THE NEXT LINE TO SEARCH *AROUND* THE OBSERVATION FREQUENCY YOU HAVE CHOSEN                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
rounding = int(np.floor(((obs_freq-1)/100)))

# PASS THE FILENAME AND FREQUENCY OF ONLY THE OBSERVATION FILES AROUND A FREQUENCY obs_freq TO A TEXT FILE                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
os.system('vap -c freq {0}/*.dzTF | egrep "( {1}..| {2}..| {3}..)" > frequency_list.txt' .format(directory,rounding,rounding+1,rounding-1))
#os.system('vap -c rcvr {0}/*.dzTF | egrep " 430" > frequency_list.txt' .format(directory))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
#os.system('vap -c freq {0}/*.x.zap | egrep "( {1}..)" > frequency_list.txt' .format(directory,rounding))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
#os.system('vap -c rcvr {0}/*.x.zap | egrep " 430" > frequency_list.txt' .format(directory))                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       

# CREATE AN ARRAY WHICH SEPARATES THE FILENAME AND THE FREQUENCY                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
y = np.genfromtxt("frequency_list.txt", dtype=[('filename','S50'),('freq','S50')])
# NUMBER OF OBSERVATIONS AT THE DESIRED FREQUENCY                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
numrows = y.shape[0]

# SORT OBSERVATIONS INTO DATE ORDER. ASSUMES THAT CHARACTERS 2-7 ARE YYMMDD OF THE OBSERVATION                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
new_y = []

for i in range(numrows):
    new_y.append(y['filename'][i])

new_y.sort(key=lambda new_y: new_y[1:7])

print new_y

# LOOP OVER EACH EPOCH (ROW) AND PERFORM PDV TO FIND OUT HOW MANY OBSERVATIONS REMAIN THAT HAVE bins PHASE BINS - ASSUMES THAT THE OBSERVATIONS HAVE ALREADY BEEN TIME AND FREQUENCY SCRUNCHED.                                                                                                                                                                                                                                                                                                                                                                                                                                    

rem_count=0

for i in range(numrows):
    epoch_name = new_y[i]
    print "(1) CURRENTLY:",epoch_name
    os.system('pdv -ZFTt {0}/{1} > temp.txt '.format(directory,epoch_name))
    stokes_line = np.genfromtxt('temp.txt', usecols=3, dtype=[('stokesI','float')], skip_header=1)
    if len(stokes_line)!=bins:
        rem_count+=1
num_left=numrows-rem_count

# LOOP OVER EACH EPOCH (ROW) AND PERFORM PDV AND READ INTO A TEXT FILE                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

stokes_list = np.zeros((num_left,bins))
b=0
removed=0

for i in range(numrows):
    print "(2) CURRENTLY:",epoch_name
    epoch_name = new_y[i]
    os.system('pdv -ZFTt {0}/{1} > temp.txt '.format(directory,epoch_name))
    stokes_line = np.genfromtxt('temp.txt', usecols=3, skip_header=1)
    if len(stokes_line) !=bins:
        removed+=1
    else:
        os.system('vap -nc "mjd, length" {0}/{1} >> mjd.txt'.format(directory,new_y[i]))
        stokes_list[b] = stokes_line
        b+=1

# PRODUCE AN ARRAY OF OBSERVATION MJDS AND LENGTHS OF THE OBSERVATIONS                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             

new_mjd = []

mjdarray = np.genfromtxt("mjd.txt", dtype=[('filename2','S40'),('mjd','f8'),('len','f8')])

for i in range(num_left):
    new_mjd.append(mjdarray['mjd'][i])
new_mjd.sort()

mjd_length = np.zeros((num_left,2))

for i in range(num_left):
    mjd_length [i,0] = mjdarray['mjd'][i]
    mjd_length [i,1] = mjdarray['len'][i]

mjd_length_sorted = mjd_length[mjd_length[:,0].argsort()]

new_length = mjd_length_sorted[:,1]
new_mjd = mjd_length_sorted[:,0]

# CREATES THE 2D MATRIX                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            

print "Number of profiles that are not",bins,"bins and therefore removed:",removed,"out of",numrows,"original profiles"
stokes_columns = np.transpose(stokes_list)
stokes_columns2 = np.zeros((bins+2,num_left))
stokes_columns2 = np.vstack([stokes_columns,new_length,new_mjd])

np.savetxt('{0}_{1}list_{2}.txt'.format(pulsar,obs_freq,bins),stokes_columns2, delimiter='\t')

# REMOVES ANY FILES THAT WERE TEMPORARILY CREATED                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  

try:
    os.remove('temp.txt')
except os.error:
    pass

try:
    os.remove('length.txt')
except os.error:
    pass

try:
    os.remove('frequency_list.txt')
except os.error:
    pass

try:
    os.remove('mjd.txt')
except os.error:
    pass

