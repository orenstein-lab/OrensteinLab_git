from strain_control.strain_client import StrainClient
'''
low level functions for controlling strain cell power suppl via strain server. Note that for these functions to work, the strain server must be running.
'''

######################
### Core Functions ###
######################

def set_strain_ps(voltage, sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    sc.set_ps(voltage)

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
        print(pos)

    return voltage

def read_strain_capacitance(sc=None):

    sc_passed = True
    if sc==None:
        sc = initialize_strain_cell_client()
        sc_passed = False

    cap = sc.get_cap()

    if sc_passed == False:
        print(pos)

    return cap

def initialize_strain_cell_client():
    return StrainClient()

def close_strain_cell_client(obj):
    return 0
