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
import time
import requests

IP = CONFIG_DICT['Orpheus IP']
PORT_T1 = CONFIG_DICT['Orpheus Port T1']
PORT_T2 = CONFIG_DICT['Orpheus Port T2']
SERIAL_T1 = CONFIG_DICT['Orpheus Serial T1']
SERIAL_T2 = CONFIG_DICT['Orpheus Serial T2']

class Orpheus:

    def __init__(self, ip, port, serial):
        self.base = (
            f"http://{ip}:{port}/{serial}/v0/PublicAPI"
        )

    def get(self, endpoint):
        r = requests.get(
            self.base + endpoint,
            timeout=5
        )
        r.raise_for_status()
        return r.json()

    def put(self, endpoint, data):
        r = requests.put(
            self.base + endpoint,
            json=data,
            timeout=5
        )
        r.raise_for_status()
        return r

    def read_wavelength(self):
        """Read the current OPA wavelength in nm."""
        url = f"{self.base}/Optical/WavelengthControl/Output/Wavelength"

        response = requests.get(url, timeout=5)
        response.raise_for_status()

        return float(response.json())


    def set_wavelength(self, wavelength_nm):

        self.put(
            "/Optical/WavelengthControl/"
            "SetWavelengthUsingAnyInteraction",
            wavelength_nm
        )

        self.wait_until_ready()

    def wait_until_ready(self):

        while True:

            state = self.get(
                "/Optical/WavelengthControl/Output"
            )

            progress = (
                100 *
                state["WavelengthSettingCompletionPart"]
            )

            print(
                f"\rSetting wavelength: "
                f"{progress:.0f}%",
                end=""
            )

            if state["IsWaitingForUserAction"]:

                print("\nManual action required:")

                for message in state["Messages"]:
                    print(message["Text"])

                raise RuntimeError(
                    "ORPHEUS waiting for user action"
                )

            if not state[
                "IsWavelengthSettingInProgress"
            ]:
                print("\nWavelength set.")
                break

            time.sleep(0.2)

    def shutter_open(self):
        self.put(
            "/ShutterInterlock/OpenCloseShutter",
            True
        )

    def shutter_close(self):
        self.put(
            "/ShutterInterlock/OpenCloseShutter",
            False
        )

    def shutter_state(self):
        return self.get(
            "/ShutterInterlock/IsShutterOpen"
        )

#
# T1
#

def read_t1(obj=None):

    obj_passed=True
    if obj==None:
        obj=init_t1()
        obj_passed=False

    ##
    ## new code goes here
    ##
    position = obj.read_wavelength()

    if obj_passed==False:
        close_t1(obj)
        return position, None
    else:
        return position, obj

def move_t1(val, obj=None, check_stability=True):

    obj_passed=True
    if obj==None:
        obj=init_t1()
        obj_passed=False

    obj.set_wavelength(val)

    if check_stability==True:
        ## check stability
        pass

    if obj_passed==False:
        close_t1(obj)
        return None
    else:
        return obj

def init_t1():
    obj = Orpheus(ip=IP, port=PORT_T1, serial=SERIAL_T1)
    return obj

def close_t1(obj):
    return 0


#
# T2
#

def read_t2(obj=None):

    obj_passed=True
    if obj==None:
        obj=init_t2()
        obj_passed=False

    ##
    ## new code goes here
    ##
    position = obj.read_wavelength()

    if obj_passed==False:
        close_t2(obj)
        return position, None
    else:
        return position, obj

def move_t2(val, obj=None, check_stability=True):

    obj_passed=True
    if obj==None:
        obj=init_t2()
        obj_passed=False

    obj.set_wavelength(val)

    if check_stability==True:
        ## check stability
        pass

    if obj_passed==False:
        close_t2(obj)
        return None
    else:
        return obj

def init_t2():
    obj = Orpheus(IP, PORT_T2, SERIAL_T2)
    return obj

def close_t2(obj):
    return 0