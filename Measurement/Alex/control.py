'''
Main file for controlling lab equipment and orchestrating measurements
'''

from strain_control.strain_client import StrainClient
import time
import numpy as np
import matplotlib.pyplot as plt

def save_data(fname, data, header, metadata=None):
    '''
    utility function for saving data to a file, with optional metadata

    args:
        - fname(string):           full path to datafile
        - data(array):             (n,m) array containing data
        - header(array):           (m) array of strings labeling each column
        - metadata(dict):          a dictionary of parameters to store at the top of the datafile under [Metadata]

    returns: None
    '''
    if not(len(header) == len(data[0,:])):
        raise ValueError('number of header items does not match number of data columns.')
    with open(fname, 'w') as f:
        if not(metadata==None):
            f.write('[METADATA]\n')
            for key, value in metadata.items():
                f.write(f'{key}:\t{value}\n')
            f.write('[DATA]\n')
        for item in header:
            f.write(str(item)+'\t')
        f.write('\n')
        for line in data:
            for item in line:
                f.write(str(item)+'\t')
            f.write('\n')

def measure_strain_cell_capacitor(fname, sc, num_points=1000, dt_min=0.1):
    '''
    mesaures the
    '''
    t = np.zeros(num_points)
    cap = np.zeros(num_points)
    t0 = time.time()
    i = 0
    t_old = t0
    while i < num_points:
        t_new = time.time()
        if t_new - t_old > dt_min:
            c = sc.get_cap()
            cap[i] = c
            t[i] = t_new - t0
            i = i + 1
            t_old = t_new
    save_data(f, np.transpose([t, cap]), ['Time', 'Cap'])
