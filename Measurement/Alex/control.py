'''
Main file for controlling lab equipment and orchestrating measurements
'''

from strain_control.strain_client import StrainClient
import time
import numpy as np
import matplotlib.pyplot as plt

def save_to_file(fname, data, header, metadata=None):
    '''
    utility function for saving data to a file, with optional metadata

    args:
        - fname(string):           full path to datafile
        - data(array):             (n,m) array containing data
        - header(array):           (m) array of strings labeling each column
        - metadata(str):           in progress

    returns: None
    '''
    if not(len(header) == len(data[0,:])):
        raise ValueError('number of header items does not match number of data columns.')
    with open(fname, 'w') as f:
        for item in header:
            f.write(str(item)+'\t')
        f.write('\n')
        for line in data:
            for item in line:
                f.write(str(item)+'\t')
            f.write('\n')


filename = 'G:/Shared drives/Orenstein Lab/Data/alex/20220818_strain_cell_capacitor_noise_25ft_coax_300kHz.txt'
num_points = 1000
t = np.zeros(num_points)
cap = np.zeros(num_points)
t0 = time.time()
i = 0
t_old = t0
while i < num_points:
    t_new = time.time()
    if t_new - t_old > 0.1:
        c = sc.get_cap()
        cap[i] = c
        t[i] = t_new - t0
        i = i + 1
        t_old = t_new
save_to_file(filename, np.transpose([t, cap]), ['Time', 'Cap'])
