from sklearn import linear_model
import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt
import numpy as np
import csv
import pandas as pd
from lib.heat_map import heatmap                # Custom library to plot param_nametmap

'''
This Script read data from a csv file, where each row define one materials and each column contain different parameters. Column #20 represents the proton migration energy barrier.
We do a linear regression analysis between different columns to find which of the parameters are directly related to energy barrier of proton migration.
At last we plot the scores for different stractural parameters for provide a visual understanding which of the parameter proton migration relates to.  
'''

reader = csv.reader(open('./data3.csv', 'rU'), dialect=csv.excel_tab) # File name, where the structure data and proton migration energy barrier for each materials are located
raw_data = []
param_nameder = []
row_count=0
for line in reader:
    if row_count == 0:
        row_count = 1
        param_nameder = line[0].split(',')
        continue
    raw_data.append (line[0].split(','))                   # Stores data from csv file


#  Convert csv data as panda frame
df_ = pd.DataFrame(np.array(raw_data),columns = param_nameder)

regra = [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19, 23, 24, 25, 26, 27, 28, 29, 35, 37, 38, 42, 51]  # Columns for analysis for linear regression
target = df_.iloc[:, [20]]
y = target  # Stores target column



ht_temp = []
ht = []
param_name = []
for column_iter_1 in regra:
    for column_iter_2 in regra:
        df = df_.iloc[:, [column_iter_1, column_iter_2]]                        # Set the column for regression analysis
        X = df
        lm = linear_model.LinearRegression()            # Define the linear model from sklearn
        model = lm.fit(X, y)                            # Linear fit the model with
        score = lm.score(X, y) ** 2                     # Calculate the regression score
        ht_temp.append(score)                           # Temporarily stores score value for one set of regression.
    ht.append(ht_temp)                                  # Store score in a array for plotting
    ht_temp = []
    param_name.append(param_nameder[column_iter_1])

fig, ax = plt.subplots()                                # Define subplot for plotting
im, cbar = heatmap.draw (np.array(ht), param_name, param_name, ax = ax, cmap = "YlGn")
fig.tight_layout()
plt.show()                                              # Plot the param_namet map
