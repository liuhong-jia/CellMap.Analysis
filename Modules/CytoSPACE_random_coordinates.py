## Generate single-cell random coordinates based on CytoSPACE results

import pandas as pd
import numpy as np
assigned_locations = pd.read_csv('assigned_locations.csv')

max_num_cells=50000
geometry="honeycomb"

if assigned_locations.shape[0] > max_num_cells:
        assigned_locations = assigned_locations.sample(max_num_cells)

def rand_jitter(arr,interval):
    return arr + np.random.uniform(-interval/2,interval/2,len(arr))

X = assigned_locations.iloc[:,-2]
Y = assigned_locations.iloc[:,-1]
cell_types = assigned_locations['CellType'].values

scale = Y.max() < 500 and ((Y - Y.round()).abs() < 1e-5).all()
y_int = 1 if scale else np.median(np.unique(np.diff(np.sort(np.unique(Y)))))
x_int = 1 if scale else np.median(np.unique(np.diff(np.sort(np.unique(X)))))

y_interval = y_int
x_interval = x_int

X = rand_jitter(X.values,x_interval)
Y = rand_jitter(Y.values,y_interval)


##将X和Y添加到assigned_locations 中

coordinates_df = pd.DataFrame({'cell.x': X, 'cell.y': Y})

assigned_locations['cell.x'] = coordinates_df['cell.x']
assigned_locations['cell.y'] = coordinates_df['cell.y']

assigned_locations.to_csv('locations.csv', index=False)