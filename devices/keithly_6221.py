'''
Keithly 6221
'''
from OrensteinLab_git.configuration import CONFIG_DICT
import serial
import time
import pickle
import os

COM_PORT = CONFIG_DICT['Keithly COM Port']
KEITHLY_HANDLE_FNAME = CONFIG_DICT['Keithly Handle']

def read_frequency(obj=None):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    obj.write(f'SOUR:WAVE:FREQ?\r'.encode())
    obj.flush()
    position = float(obj.readline().decode('ascii').strip())

    if obj_passed==False:
        close(obj)
        return position, None
    else:
        return position, obj

def set_frequency(val, obj=None, check_stability=True, wait_time=0):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    disable_output(obj)
    obj.write(f'SOUR:WAVE:FREQ {val}\r'.encode())
    obj.flush()
    enable_output(obj)
    time.sleep(3/val)
    time.sleep(wait_time)


    if check_stability==True:
        ## check stability
        pass

    if obj_passed==False:
        close(obj)
        return None
    else:
        return obj

def read_amplitude(obj=None):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    obj.write(f'SOUR:WAVE:AMPL?\r'.encode())
    obj.flush()
    position = float(obj.readline().decode('ascii').strip())

    if obj_passed==False:
        close(obj)
        return position, None
    else:
        return position, obj

def set_amplitude(val, obj=None, check_stability=True, best_range=True, wait_time=1):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    disable_output(obj)
    obj.write(f'SOUR:WAVE:AMPL {val}\r'.encode())
    if best_range:
        obj.write(f'SOUR:WAVE:RANG BEST\r'.encode())
    obj.flush()
    enable_output(obj)
    time.sleep(wait_time)


    if check_stability==True:
        ## check stability
        pass

    if obj_passed==False:
        close(obj)
        return None
    else:
        return obj

def read_range(obj=None):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    obj.write(f'SOUR:CURR:RANG?\r'.encode())
    obj.flush()
    position = float(obj.readline().decode('ascii').strip())

    if obj_passed==False:
        close(obj)
        return position, None
    else:
        return position, obj

def set_range(val, obj=None, check_stability=True):
    # range in amps

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    obj.write(f'SOUR:CURR:RANG {val}\r'.encode())
    obj.flush()

    if check_stability==True:
        ## check stability
        pass

    if obj_passed==False:
        close(obj)
        return None
    else:
        return obj

def read_voltage_compliance(obj=None):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    obj.write(f'CURR:COMP?\r'.encode())
    obj.flush()
    position = float(obj.readline().decode('ascii').strip())

    if obj_passed==False:
        close(obj)
        return position, None
    else:
        return position, obj

def set_voltage_compliance(val, obj=None, check_stability=True):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    obj.write(f'CURR:COMP {val}\r'.encode())
    obj.flush()

    if check_stability==True:
        ## check stability
        pass

    if obj_passed==False:
        close(obj)
        return None
    else:
        return obj

def init():
    obj = SerialManager.get_connection()
    return obj

def close(obj):
    obj.close()

def enable_output(obj=None):
    # enables Keithly output

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    obj.write(b'SOUR:WAVE:ARM\r')
    obj.flush()
    obj.write(b'OUTP ON\r')
    obj.flush()
    obj.write(b'SOUR:WAVE:INIT\r')
    obj.flush()

    if obj_passed==False:
        close(obj)
        return None
    else:
        return obj

def disable_output(obj=None):
    # disables Keithly output
    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    obj.write(b'OUTP OFF\r')

    if obj_passed==False:
        close(obj)
        return None
    else:
        return obj


class SerialManager:
    _connection = None

    @classmethod
    def get_connection(cls):
        # Reuse existing open connection if it matches
        if cls._connection is not None and cls._connection.is_open:
            return cls._connection

        cls._connection = serial.Serial(port=COM_PORT,             # Replace with your serial port
                    baudrate=9600,           # Set the baud rate for your device
                    parity=serial.PARITY_NONE,
                    stopbits=serial.STOPBITS_ONE,  # 1 stop bits
                    bytesize=serial.EIGHTBITS,
                    timeout=1,                # Timeout for read operations
                    xonxoff=False,
                    rtscts=False
                    )
        return cls._connection

    @classmethod
    def close(cls):
        if cls._connection and cls._connection.is_open:
            cls._connection.close()
        cls._connection = None
