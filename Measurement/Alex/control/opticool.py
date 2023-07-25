import sys
import telnetlib
from OrensteinLab_git.configuration import config_dict

host = config_dict['Opticool IP']
port = config_dict['Opticool Port']
OPTICOOL_HANDLE_NAME = config_dict['Opticool Handle']

######################
### Core Functions ###
#####################

def set_opticool_temperature(temperature, optc=None, rate=10, mode=0):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

        message=('TEMP '+str(set_point)+', '+str(rate)+', '+str(mode)).encode('ascii')
        optc.write(message+('\r\n').encode('ascii'))
        output=optc.read_some().decode('ascii')

    if obj_passed==False:
        close_opticool(optc)

    return output

def read_opticool_temperature(optc=None):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    message=('TEMP?').encode('ascii')
    optc.write(message+('\r\n').encode('ascii'))
    output=optc.read_some().decode('ascii')
    output_value = output.split(',')[1].strip()
    output_status = output.split(',')[3].strip()
    output_status = output_status.replace("\"","")
    if output_value == 'nan':
        temp = output_value
    else:
        temp = float(output_value)

    if obj_passed==False:
        close_opticool(optc)

    return temp

def set_opticool_field(field, optc=None, rate=110, mode=0):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    message=('FIELD '+str(set_point)+', '+str(rate)+', '+str(mode)+', 1').encode('ascii')
    optc.write(message+('\r\n').encode('ascii'))
    output=optc.read_some().decode('ascii')

    if obj_passed==False:
        close_opticool(optc)

    return output

def read_opticool_field(optc=None):

    obj_passed=True
    if optc==None:
        optc = initialize_opticool()
        obj_passed==False

    message=('FIELD?').encode('ascii')
    optc.write(message+('\r\n').encode('ascii'))
    output=optc.read_some().decode('ascii')
    output_value = output.split(',')[1].strip()
    output_status = output.split(',')[3].strip()
    output_status = output_status.replace("\"","")
    if output_value == 'nan':
        field = output_value
    else:
        field = float(output_value)

    if obj_passed==False:
        close_opticool(optc)

    return field

def initialize_opticool():
    try:
        with open(OPTICOOL_HANDLE_FNAME, 'rb') as f:
            optc = pickle.load(f)
    except: # if there is not already a connection, above fails and we instantiate a new handle
        optc=telnetlib.Telnet(host, port, timeout=15)
        optc.read_until(('Connected to QDInstrument Socket Server.\r\n').encode('ascii'))
        with open(OPTICOOL_HANDLE_FNAME, 'wb') as f:
            pickle.dump(optc, f)
    return optc

def close_opticool(optc):
    optc.close()
