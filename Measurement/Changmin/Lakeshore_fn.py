import OrensteinLab.Instrument.Lakeshore.Lakeshore335 as ls
import time
from IPython.display import clear_output

lsobj = ls.initialization_lakeshore335()

def read_temperature():
    return ls.read_temperature(lsobj)

def set_ramp(ramprate):
    ls.set_ramp(lsobj, 1, 1, ramprate)

def set_setpoint(Tset):
    ls.set_setpoint(lsobj, 1, Tset)

def stab_temperature(Tset,Ttol,stabtime):
    set_ramp(0)
    time.sleep(0.1)
    set_setpoint(Tset)
    time.sleep(0.1)
    t0 = time.time()
    tnow = time.time()
    thisT = float(read_temperature())
    time.sleep(0.1)
    while tnow - t0 < stabtime or abs(thisT - Tset) > Ttol:
        clear_output(wait=True)
        print(thisT)
        if abs(thisT - Tset) > Ttol:
            t0 = time.time()
        time.sleep(0.1)
        tnow = time.time()
        thisT = float(read_temperature())
    print('Temperature stabilized!')
