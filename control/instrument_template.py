'''
instrument_template -  each instrument at minimum must contain

    - init(**kwargs) - returns instrument handle object
    - read(obj=None, **kwargs) - return dictionary {key:value, ...} where keys are name of measured quantity and values are values of that quantity
    - close(obj, **kwargs) - properly closes instrument handle given object

    Other requirements:
    - read functions must auto-initialize and close if no instrument handle is passed.

beyond this, users may write any code for reading instruments. A given file may also contain different instruments associated with a single device, and using wrappers for repeated actions is highly encouraged.
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

def init():
    return 0

def close(obj):
    return 0