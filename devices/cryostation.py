'''
file for Montana cryostation
'''
from OrensteinLab_git.configuration import CONFIG_DICT
import OrensteinLab_git.devices.montana_cryostation.Commands as cmd
import OrensteinLab_git.devices.montana_cryostation.CryostationComm as com
import time


IP = CONFIG_DICT['Montana Cryostation IP']
PORT = CONFIG_DICT['Montana Cryostation Port']

def read_cryo_field(obj=None):

    obj_passed=True
    if obj==None:
        obj=init_cryo()
        obj_passed=False

    position = float(cmd.get_magnet_target_field(obj)[1])

    if obj_passed==False:
        close_cryo(obj)
        return position, None
    else:
        return position, obj

def move_cryo_field(val, obj=None, check_stability=True, wait_time=0):

    obj_passed=True
    if obj==None:
        obj=init_cryo()
        obj_passed=False

    cmd.set_target_field(obj, val)

    if check_stability==True:
        ## check stability
        pass

    time.sleep(wait_time)

    if obj_passed==False:
        close_cryo(obj)
        return None
    else:
        return obj

def init_cryo():
    obj = com.CryoComm(IP, PORT)
    return obj

def close_cryo(obj):
    # if obj.socket:
    #     obj.socket.shutdown(1)
    #     obj.socket.close()
    del obj
    #obj.__del__()
    return 0
