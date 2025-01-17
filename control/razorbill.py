from strain_control.strain_client import StrainClient
import numpy as np
import time
'''
low level functions for controlling strain cell power suppl via strain server. Note that for these functions to work, the strain server must be running.
'''

######################
### Core Functions ###
######################

def set_strain_ps(voltage, sc=None, wait_time=0, check_stability=True, tol=0.1):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    sc.set_ps(voltage)
    if check_stability:
        while np.abs(read_strain_ps(sc)-voltage) >= tol:
            time.sleep(0.1)
    time.sleep(wait_time)

def set_voltage_1(voltage, sc=None, wait_time=0, check_stability=True, tol=0.1):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    sc.set_voltage(1, voltage)
    if check_stability:
        while np.abs(read_voltage_1(sc)-voltage) >= tol:
            time.sleep(0.1)
    time.sleep(wait_time)

def set_voltage_2(voltage, sc=None, wait_time=0, check_stability=True, tol=0.1):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    sc.set_voltage(2, voltage)
    if check_stability:
        while np.abs(read_voltage_2(sc)-voltage) >= tol:
            time.sleep(0.1)
    time.sleep(wait_time)

def set_strain_capacitance(cap, sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    sc.set_cap(cap)

def read_strain_ps(sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    voltage = sc.get_ps()

    if sc_passed == False:
        print(voltage)

    return voltage

def read_voltage_1(sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    voltage = sc.get_voltage(1)

    if sc_passed == False:
        print(voltage)

    return voltage

def read_voltage_2(sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    voltage = sc.get_voltage(2)

    if sc_passed == False:
        print(voltage)

    return voltage

def read_strain_capacitance(sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    cap = sc.get_cap()

    if sc_passed == False:
        print(cap)

    return cap

def initialize_strain_cell_client():
    return StrainClient()

def close_strain_cell_client(obj):
    return 0