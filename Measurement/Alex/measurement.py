'''
Main file for controlling lab equipment and orchestrating measurements, with a specific eye to procedures
'''

from strain_control.strain_client import StrainClient
import OrensteinLab_git.Measurement.Alex.control as ctrl
import OrensteinLab_git.Instrument.montana.cryocore as cryocore
import time
import numpy as np
import matplotlib.pyplot as plt

def save_data_to_file(fname, data, header, metadata=None):
    '''
    utility function for saving data to a file, with optional metadata

    args:
        - fname(string):           full path to datafile
        - data(array):             (n,m) array containing data
        - header(array):           (m) array of strings labeling each column
        - metadata(dict):          a dictionary of parameters to store at the top of the datafile under [Metadata]

    returns: None
    '''
    if not(len(header) == len(data[0,:])):
        raise ValueError('number of header items does not match number of data columns.')
    with open(fname, 'w') as f:
        if not(metadata==None):
            f.write('[METADATA]\n')
            for key, value in metadata.items():
                f.write(f'{key}:\t{value}\n')
            f.write('[DATA]\n')
        for item in header:
            f.write(str(item)+'\t')
        f.write('\n')
        for line in data:
            for item in line:
                f.write(str(item)+'\t')
            f.write('\n')

def measure_strain_cell_capacitor(fname, sc, num_points=1000, dt_min=0.1):
    '''
    mesaures the
    '''
    t = np.zeros(num_points)
    cap = np.zeros(num_points)
    t0 = time.time()
    i = 0
    t_old = t0
    while i < num_points:
        t_new = time.time()
        if t_new - t_old > dt_min:
            c = sc.get_cap()
            cap[i] = c
            t[i] = t_new - t0
            i = i + 1
            t_old = t_new
    save_data_to_file(fname, np.transpose([t, cap]), ['Time', 'Cap'])

def strain_cell_temperature_calibration(fname1, fname2, filename_head, sc, cryo, temps, lakeshore_stability, cap_stability, cryo_stability_high, cryo_stability_low, mode, wait_time=1):
    '''
    runs a cooldown and warmup of Montana CryoAdvance and monitors platform temperature vs strain cell capacitance. The strain cell should be loaded with a titanium dummy sample. For consistency, block temperature should be measured NOT sample temperature, which will change from setup to setup.

    args:
        - fname1:           file path to a file to log continuosly during a ramp
        - fname2:           file path to a file to log only at times when both temperature and capacitor are stable
        - sc:               StrainClient object
        - cryo:             Montana CryoCore object
        - temps:            list of temperatures. Note that Montana seems to like integer values.
        - mode:             1 for keeping cryostat at base and using heater and 2 for changing cold plate temp
    '''
    filename1 = filename_head+'\\'+fname1+'.dat'
    filename2 = filename_head+'\\'+fname2+'.dat'
    with open(filename1, 'a') as f1:
        f1.write('Time' + '\t' + 'Setpoint Temperature (K)' + '\t' + 'Platform Temperature (K)' + '\t' + 'Lakeshore Temperature (K)' + '\t' + 'Capacitance' + '\n')
    with open(filename2, 'a') as f2:
        f2.write('Platform Temperature (K)' + '\t' + 'Lakeshore Temperature (K)' + '\t' + 'Capacitance' + '\n')

    t0 = time.time()
    cryo.cooldown()
    print('Cooling down cryostat.')

    if mode==1: # cooldown to base temperature and then set temperature with lakeshore, waiting for both lakeshore and capacitance to stabilize before moving on.
        print('Measuring in mode 1.')
        setpoints = np.sort(temps) # measure from base to high temp.
        target_stability = cryo_stability_low
        cryo.set_platform_target_temperature(4.63)
        cryo.set_platform_stability_target(target_stability)
        ctrl.set_temperature(0)
        while True:
            time.sleep(wait_time)
            stability_ok, is_stable = cryo.get_platform_temperature_stable()
            if is_stable:
                print(f'Stabilized platform to base temperature.')
                break
        cryo.set_platform_target_temperature(0)
        for sp in setpoints:
            print(f'Setting setpoint to {sp} K')
            ctrl.set_temperature(sp)
            lakeshore_temps = []
            while True:
                time.sleep(wait_time)
                lakeshore_temps.append(ctrl.read_temperature())
                with open(filename1, 'a') as f1:
                    f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                         str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                         str(format(float(lakeshore_temps[-1]), '.5f')) + '\t' +
                         str(format(float(sc.get_cap()), '.5f')) + '\n')
                if len(lakeshore_temps) > 120:
                    mean = np.mean(np.asarray(lakeshore_temps[-120:]))
                    std = np.std(np.asarray(lakeshore_temps[-120:]))
                    if ((std < lakeshore_stability) or (len(lakeshore_temps) > 7200)) and ((mean > sp-0.01) and (mean < sp+0.01)):
                        if std < lakeshore_stability:
                            print(f'Stabilized Lakeshore temperature at {mean} K')
                        if len(lakeshore_temps) > 7200:
                            print('Exceeded maximum soak time')
                        print(f'Lakeshore noise: {std}')
                        caps = []
                        while True:
                            time.sleep(wait_time)
                            caps.append(sc.get_cap())
                            lakeshore_temps.append(ctrl.read_temperature())
                            with open(filename1, 'a') as f1:
                                f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                                     str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                     str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                     str(format(float(caps[-1]), '.5f')) + '\n')
                            if len(caps)>300: # mininum soak time of 5 minutes
                                if (np.std(np.asarray(caps[-300:])) < cap_stability) or (len(caps) > 7200):
                                    print('STD: '+str(np.std(np.asarray(caps[-300:]))))
                                    if (np.std(np.asarray(caps[-300:])) < cap_stability):
                                        print('Stdev below accepted value')
                                    elif len(caps) > 7200: # maximum soak time of 2 hours
                                        print('Exceeded maximum soak time')
                                    print(f'Stabilized capacitance measurement, writing to file')
                                    with open(filename2, 'a') as f2:
                                        f2.write(str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                             str(format(np.mean(lakeshore_temps[-15:]), '.5f')) + '\t' +
                                             str(format(np.mean(caps[-15:]), '.5f')) + '\n')
                                    break
                                elif (np.std(np.asarray(caps[-300:])) < cap_stability*1.1) and (len(caps) > 3600):
                                    print('STD: '+str(np.std(np.asarray(caps[-300:]))))
                                    print('Stdev within 10 percent of accepted value after 1 hour')
                                    print(f'Stabilized capacitance measurement, writing to file')
                                    with open(filename2, 'a') as f2:
                                        f2.write(str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                             str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                             str(format(np.mean(caps[-15:]), '.5f')) + '\n')
                                    break
                        break

    if mode==2: # on cooldown, set cryostat target temperature and wait for lakeshore and capacitance - this still needs some work to work around stability issues at 10-12K by switching the stability criteria for both cryostat and lakeshore.
        print('Measuring in mode 1.')
        setpoint = np.flip(np.sort(temps)) # measure from high temp to low temp.
        for sp in setpoints:
            print(f'Setting setpoint to {sp} K')
            if sp >= 10:
                target_stability = cryo_stability_high
            else:
                target_stability = cryo_stability_low
            cryo.set_platform_target_temperature(int(sp))
            cryo.set_platform_stability_target(target_stability)
            while True:
                time.sleep(wait_time)
                with open(filename1, 'a') as f1:
                    f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                         str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                         str(format(float(ctrl.read_temperature()), '.5f')) + '\t' + str(format(float(sc.get_cap()), '.5f')) + '\n')
                stability_ok, is_stable = cryo.get_platform_temperature_stable()
                if is_stable:
                    print(f'Stabilized platform temperature at {cryo.get_platform_temperature()[1]} K')
                    lakeshore_temps = []
                    while True:
                        time.sleep(wait_time)
                        lakeshore_temps.append(ctrl.read_temperature())
                        with open(filename1, 'a') as f1:
                            f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                                 str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                 str(format(float(lakeshore_temps[-1]), '.5f')) + '\t' +
                                 str(format(float(sc.get_cap()), '.5f')) + '\n')
                        if len(lakeshore_temps) > 120:
                            if (np.std(np.asarray(lakeshore_temps[-120:])) < lakeshore_stability) or (len(lakeshore_temps) > 7200):
                                print(f'Stabilized Lakeshore temperature at {ctrl.read_temperature()} K')
                                print(f'Lakeshore noise: {np.std(np.asarray(lakeshore_temps[-120:]))}')
                                caps = []
                                while True:
                                    time.sleep(wait_time)
                                    caps.append(sc.get_cap())
                                    with open(filename1, 'a') as f1:
                                        f1.write(str(format(float(time.time()-t0), '.5f')) + '\t' + str(sp) + '\t' +
                                             str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                             str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                             str(format(float(caps[-1]), '.5f')) + '\n')
                                    if len(caps)>300: # mininum soak time of 5 minutes
                                        if (np.std(np.asarray(caps[-300:])) < cap_stability) or (len(caps) > 7200):
                                            print('STD: '+str(np.std(np.asarray(caps[-300:]))))
                                            if (np.std(np.asarray(caps[-300:])) < cap_stability):
                                                print('Stdev below accepted value')
                                            elif len(caps) > 7200: # maximum soak time of 2 hours
                                                print('Exceeded maximum soak time')
                                            print(f'Stabilized capacitance measurement, writing to file')
                                            with open(filename2, 'a') as f2:
                                                f2.write(str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                                     str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                                     str(format(np.mean(caps[-15:]), '.5f')) + '\n')
                                            break
                                        elif (np.std(np.asarray(caps[-300:])) < cap_stability*1.1) and (len(caps) > 3600):
                                            print('STD: '+str(np.std(np.asarray(caps[-300:]))))
                                            print('Stdev within 10 percent of accepted value after 1 hour')
                                            print(f'Stabilized capacitance measurement, writing to file')
                                            with open(filename2, 'a') as f2:
                                                f2.write(str(format(float(cryo.get_platform_temperature()[1]), '.5f')) + '\t' +
                                                     str(format(float(ctrl.read_temperature()), '.5f')) + '\t' +
                                                     str(format(np.mean(caps[-15:]), '.5f')) + '\n')
                                            break
                                break
                    break

    cryo.warmup()
    print('Warming up cryostat')

def zero_strain_cell(sc, slew_rate=1, target_voltage=120, tol=0.1):
    '''
    carries out strain cell zeroing procedure as laid out in Razorbill documentatin. Energise both inner and outer (channel 1 and 2) stacks to +120V at room temperature with NO SAMPLE MOUNTED, and then allow stacks to return slowly to 0V by setting votlage to 0V and slew rate to 0.1, then turning off the outputs.

    args:
        - sc:           StrainClient object

    returns:
        - cap:          0 strain measurement of capacitance.
    '''
    sc.set_slew_rate(slew_rate)
    sc.set_output(1,1)
    sc.set_output(2,1)
    sc.set_voltage(1, target_voltage)
    sc.set_voltage(2, target_voltage)
    cond = True
    while cond:
        time.sleep(0.1)
        v1 = sc.get_voltage(1)
        v2 = sc.get_voltage(2)
        if (v1 > target_voltage-tol and v1 < target_voltage+tol) and (v2 > target_voltage-tol and v2 < target_voltage+tol):
            print('Reached 120 V on both channels')
            cond = False
    time.sleep(4)
    sc.set_slew_rate(0.1)
    sc.set_voltage(1,0)
    sc.set_voltage(2,0)
    cond = True
    while cond:
        time.sleep(0.1)
        v1 = sc.get_voltage(1)
        v2 = sc.get_voltage(2)
        if (v1 > -tol and v1 < tol) and (v2 > -tol and v2 < tol):
            cond=False
    sc.set_output(1,0)
    sc.set_output(2,0)
    cap = 0
    for i in range(100):
        cap = cap + sc.get_cap()
        time.sleep(0.1)
    return cap/100
