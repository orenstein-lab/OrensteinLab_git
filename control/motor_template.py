'''
motor_template -  each motor at minimum must contain

    - init(**kwargs) - returns motor handle object
    - read(obj=None, **kwargs) - return motor value
    - move(val, obj=None, check_stability=True, **kwargs) - moves motor to val
    - close(obj, **kwargs) 

    Other requirements:
    - read and move functions must auto-initialize and close if no motor handle is passed.
    - move functions must check stability of motor if check_stability=True

beyond this, users may write any code for motor control. A given file may also contain different motors associated with a single device, and using wrappers for multiaxis motors which are controlled in equivalent manner is highly encouraged. all of these function can be renamed.
'''

def read(obj=None):

    obj_passed=True
    if obj==None:
        obj=init()
        obj_passed=False

    ##
    ## new code goes here    
    ##

    if obj_passed==False:
        close(obj)

    return 0

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

    return 0

def init():
    return 0

def close(obj):
    return 0
