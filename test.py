import os, glob

for filename in glob.glob("buffer*" or "OSMroads*"):
    os.remove(filename) 

