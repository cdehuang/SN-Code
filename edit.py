import os
from os.path import join, getsize
import re
import fileinput

filenames = []
#originally this walked between the directories by itself, but it would fail somewhere and was taking a long time to fix, so I just got lazy and changed the directory names.
#put a spectra.out.6.1 file in the research/Sn project directory to play with don't know why it's not working, just decided to not make it walk from lanl.ccn.spectra cause it's going wrong somewhere
for root, dirs, files in os.walk('Lanl.ccsn.spectra/z40G'):
    #gets the path of the files
    for name in files:
        filenames.append(os.path.join(root, name))
        #print filenames
for line in fileinput.input(filenames, inplace=1):
    line = re.sub(r'(?<!E)-(?=\d{3})', 'E-', line)
    print line,
