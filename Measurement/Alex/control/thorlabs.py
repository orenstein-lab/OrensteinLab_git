'''
Thorlabs
'''
from pylablib.devices import Thorlabs
from OrensteinLab_git.configuration import CONFIG_DICT
MFF_PORT = CONFIG_DICT['MFF_PORT']

def read(obj=None):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##
    position = 0

    if obj_passed==False:
        close(obj)
        return 0
    else:
        return position, obj

def move(val, obj=None, check_stability=True):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here
    ##

    if check_stability==True:
        ## check stability
        pass

    if obj_passed==False:
        close(obj)
    else:
        return obj

def init():
    return Thorlabs.MFF101(MFF_PORT)

def close(obj):
    return 0
