'''
motor_template -  each motor at minimum must contain following functions with given form of arguments

    - init(**kwargs) - returns motor handle object
    - read(obj=None, **kwargs) - return motor value and obj if it has been passed
    - move(val, obj=None, check_stability=True, **kwargs) - moves motor to val and return obj if it has been passed
    - close(obj, **kwargs)

    Other requirements:
    - read and move functions must auto-initialize and close if no motor handle is passed.
    - move functions must check stability of motor if check_stability=True

beyond this, users may write any code for motor control. A given file may also contain different motors associated with a single device, and using wrappers for multiaxis motors which are controlled in equivalent manner is highly encouraged. all of these function can be renamed.
'''
from OrensteinLab_git.configuration import CONFIG_DICT

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
        return position, None
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
        return None
    else:
        return obj

def init():
    return 0

def close(obj):
    return 0
