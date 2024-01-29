from OrensteinLab_git.Measurement.Alex.control.montana import cryocore
from OrensteinLab_git.configuration import config_dict
import time
import numpy as np

address = config_dict['Montana Address']

def set_montana_temperature(temperature, cryo=None, check_stability=True):
    '''
    sets lakeshore setpoint, waits until temperature is within tolerance of setpoint, and waits for soak time before returning.

    args:
        - temperature(float)

    returns: None
    '''
    # Lakeshore initialization
    cryo_passed = True
    if cryo is None:
        cryo = initialize_montana()
        cryo_passed = False

    temp=float(temperature)
    cryo.set_platform_target_temperature(temp)
    time.sleep(0.1)

    # check for stability
    if check_stability==True:
        is_stable = False
        while not is_stable:
            stability_ok, is_stable = cryo.get_platform_temperature_stable()
            time.sleep(5)

    if cryo_passed == False:
        close_montana(cryo)

def read_montana_temperature(cryo=None):
    '''
    reads temperature from lakeshore controller

    args: None

    returns:
        - temp(float):  read temperature
    '''
    cryo_passed = True
    if cryo is None:
        cryo = initialize_montana()
        cryo_passed = False

    temp = cryo.get_platform_temperature()[1]

    if cryo_passed == False:
        close_montana(cryo)

    return temp

def initialize_montana():
    return cryocore.CryoCore(address)

def close_montana(cryo):
    return 0
