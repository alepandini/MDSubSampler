from sklearn import preprocessing
import numpy as np
import pandas as pd
import csv
import re
import configparser

config = configparser.ConfigParser()
config.read("SAMPLE1.INI",)
path = config['FILEPATH']


trajectory_filename = path["BasePath"]+path["trajectory"]
topology_filename = path["BasePath"]+path["topology"]
   


Cov = pd.read_table(topology_filename, 
                    header=0, 
                    skipfooter=1,
                    delim_whitespace=True,
                    names=["A", "B", "C", "D", "E", "F", "G","H","I"]
                    )
print(Cov)
Dropresult = Cov.drop(['A', 'B'], axis=1)

print(Dropresult)
d = preprocessing.normalize(Dropresult, axis=0)
scaled_df = pd.DataFrame(d, columns=["C", "D", "E", "F", "G","H","I"])

scaled_df.head()