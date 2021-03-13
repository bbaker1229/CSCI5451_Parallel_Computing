#!/usr/bin/env python3

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans

data = pd.read_csv("/export/scratch/users/csci5451/pollution_Vsmall.csv")
data = data.drop(axis=1, columns=' ')

kmean = KMeans(n_clusters=5)
kmean.fit(data)
# print(kmean.cluster_centers_)
vals = np.unique(kmean.labels_, return_counts=True)
for i in vals[0]:
    print(str(i) + ": " + str(vals[1][i]))
