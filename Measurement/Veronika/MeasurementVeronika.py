import os
import sys
import numpy as np
import time
import lakeshore
import zhinst.utils as ziutils
import numpy as np
import time
from IPython.display import clear_output
import matplotlib.pyplot as plt
from IPython.display import display
def add_unique_postfix(fn):
    path, name = os.path.split(fn)
    name, ext = os.path.splitext(name)

    make_fn = lambda i: os.path.join(path, '%s_%04d%s' % (name, i, ext))

    for i in range(1, 1000):
        uni_fn = make_fn(i)
        if not os.path.exists(uni_fn):
            return uni_fn

    return None



device_id = 'DEV3537'
apilevel = 6
(daq, device, props) = ziutils.create_api_session(device_id, apilevel)

import instruments.newport as newport
controller = newport.NewportESP301.open_serial(port="COM8", baud=921600)


def initialization_lakeshore335(baud):
    return lakeshore.model_335.Model335(baud)
def read_temperature(inst):
    return float(inst.query('KRDG?'))
def read_setpoint(inst):
    return float(inst.query('SETP?'))

def read_setpoint(inst, output):
    return float(inst.query('SETP? '+str(output)))

def read_range(inst, output):
    return float(inst.query('RANGE? '+str(output)))

def set_range(inst, output, rang):
    inst.command('RANGE '+str(output)+','+ str(float(rang)))

def set_setpoint(inst, output, set_temperature):
    inst.command('SETP '+str(output)+','+str(float(set_temperature)))

def read_ramp(inst):
    output = inst.query('RAMP?').split(',')
    on_off = bool(int(output[0]))
    rate = float(output[1])
    return [on_off, rate]

def set_ramp(inst, output, on_off, rate):
    inst.command("RAMP "+str(output)+','+str(int(on_off))+','+str(rate))


def close_lakeshore335(inst):
    inst.disconnect_usb()
    


axis2_rot = newport.NewportESP301Axis(controller,2-1)
axis2_rot.enable()

axis1_rot = newport.NewportESP301Axis(controller,1-1)
axis1_rot.enable()

#pip install pyanc350
path1=os.curdir
os.chdir('C:\\Users\\Admin\\Google Drive\\Python Programs')
from pyanc350.v4 import Positioner
os.chdir(path1)
ax = {'x':0,'y':1,'z':2}
#define a dict of axes to make things simpler

anc = Positioner()

#anc.connect()

def read_attocube(axis):
    time.sleep(0.1)
    return anc.getPosition(ax[axis])*1000000



def move_attocube(inst, axis,position, tolerance, printout, stop): 
    position=position*1e-6;
    tolerance=tolerance*1e-6;
    inst.setTargetRange(ax[axis],tolerance)
    inst.setTargetPosition(ax[axis],position)
    inst.startAutoMove(ax[axis], 1, 0)

    #check what's happening
    time.sleep(0.5)
    moving = 1
    target = 0
    eotFwd=0
    eotBwd=0

    while target == 0:
        connected, enabled, moving, target, eotFwd, eotBwd, error = inst.getAxisStatus(ax[axis]) #find bitmask of status
        if eotFwd ==1:
            print('End of travel')
            break
        elif eotBwd ==1:
            print('End of travel')
            break
        elif target == 0:
            if printout:
                clear_output(wait=True)
                print('axis ' + axis + ' moving, currently at',inst.getPosition(ax[axis]))
        elif target == 1:
            if printout:
                print('axis ' + axis + ' arrived at',inst.getPosition(ax[axis]))
            if stop:
                inst.startAutoMove(ax[axis], 0, 0)

        time.sleep(0.1)


 #Plot
import matplotlib.colors as colors
def plotMap(x_start, x_end, y_start, y_end, step, demod_x, demod_y, demod_r):
    
  
    x_coordinates = np.arange(x_start, x_end+1, step)
    y_coordinates = np.arange(y_start, y_end+1, step)
    X_coor, Y_coor = np.meshgrid(x_coordinates, y_coordinates)
    
    x_num=(len(x_coordinates))
    y_num=(len(y_coordinates))
    
    demod_x0 = np.zeros((y_num, x_num))
    demod_y0 = np.zeros((y_num, x_num))
    demod_r0 = np.zeros((y_num, x_num))
    
    extent=[x_start, x_end, y_start, y_end]

    length = len(demod_x)
    y_num0 = length//x_num
    x_num0 = length-y_num0*x_num
    demod_x0[:y_num0, :] = np.reshape(demod_x[:y_num0*x_num], (y_num0, x_num))
    demod_y0[:y_num0, :] = np.reshape(demod_y[:y_num0*x_num], (y_num0, x_num))
    demod_r0[:y_num0, :] = np.reshape(demod_r[:y_num0*x_num], (y_num0, x_num))
    if (y_num0 < y_num):
        demod_x0[y_num0, :x_num0] = demod_x[y_num0*x_num:length]
        demod_y0[y_num0, :x_num0] = demod_y[y_num0*x_num:length]
        demod_r0[y_num0, :x_num0] = demod_r[y_num0*x_num:length]

    fig = plt.figure(figsize=(12,36))
    gs = fig.add_gridspec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])
    ax1.set_xlabel('x (um)')
    ax1.set_ylabel('y (um)')
   # ax1.set_xticks(x_plot_range)
   # ax1.set_yticks(y_plot_range)
    ax2.set_xlabel('x (um)')
    ax2.set_ylabel('y (um)')
    #ax2.set_xticks(x_plot_range)
    #ax2.set_yticks(y_plot_range)
    ax3.set_xlabel('x (um)')
    ax3.set_ylabel('y (um)')
    #ax3.set_xticks(x_plot_range)
    #ax3.set_yticks(y_plot_range)
    pos1 = ax1.imshow(demod_x0,cmap="Greys",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos1, ax=ax1)
    pos2 = ax2.imshow(demod_y0,cmap="Greys",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos2, ax=ax2)
    pos3 = ax3.imshow(demod_r0,cmap="Greys",origin='lower',extent=extent,norm = colors.TwoSlopeNorm(0))
    fig.colorbar(pos3, ax=ax3)
    plt.show()
    
    
def mapper (inst,xmin, xmax, ymin, ymax, step,tolerance,waittime, filename):
    filename=add_unique_postfix(filename+'.dat')
    t0=time.time()
    
    file = open(filename,'a')
    file.write('Time(s)'+'\t'+'X'+'\t'+'Y'+'\t'+'Z'+'\t'+'Rx'+'\t'+'Ry'+'Ry'+'\n')
    file.close()
    timeArray = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    xArray=np.array([])
    yArray=np.array([])
    RArray=np.array([])
    
    #take the backlash
    move_attocube(inst, 'x',xmin-10, tolerance, 1, 0)
    move_attocube(inst, 'y',ymin-10, tolerance, 1, 0)
    t0=time.time()
    #start the map
    Nx=int((xmax-xmin)/step+1)
    Ny=int((ymax-ymin)/step+1)
    xs=np.linspace(xmin, xmax, Nx)
    ys=np.linspace(ymin, ymax, Ny)
    for i in range(0, len(xs)):
        x=xs[i]
        for j in range(0, len(ys)):
            y=ys[i]
            clear_output(wait=True)
            #move
            move_attocube(inst, 'x',x, tolerance, 0, 0)
            move_attocube(inst, 'y',y, tolerance, 0, 0)
            time.sleep(waittime)
            
            #read voltage 
            sample = daq.getSample('/%s/demods/0/sample' % device)
            sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
            Rx = sample["x"][0]
            Ry = sample["y"][0]
            R = sample["R"][0]
            
            #add to arrays
            xArray = np.append(xArray, Rx)
            yArray = np.append(yArray, Ry)
            RArray = np.append(RArray, R)
            
            #time
            ts=time.time()-t0
            
            #write to file
            file = open(filename,'a')
            file.write(format(ts, '.15f')+"\t"+format(read_attocube(anc, 'x'), '.15f')+"\t"+format(read_attocube(anc, 'y'), '.15f')+"\t"+format(read_attocube(anc, 'z'), '.15f')+"\t"+format(Rx, '.15f')+'\t'+format(Ry, '.15f')+'\t'+format(R, '.15f')+'\n')
            file.close()
            #plot maps
           # plotMap(xmin, xmax, ymin, ymax, step, xArray, yArray, RArray)
            print (x, y)


def mapperZ (inst,z, xmin, xmax, ymin, ymax, step,tolerance,waittime, filename):
    filename=add_unique_postfix(filename+'.dat')
    t0=time.time()
    
    file = open(filename,'a')
    file.write('Time(s)'+'\t'+'X'+'\t'+'Y'+'\t'+'Z'+'\t'+'Rx'+'\t'+'Ry'+'\t'+'R'+'\t'+'Temp'+'\n')
    file.close()
    timeArray = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    xArray=np.array([])
    yArray=np.array([])
    RArray=np.array([])
    Temp=np.array([])
    
    #take the backlash
    move_attocube(inst, 'z',z-10, tolerance, 1, 0)
    move_attocube(inst, 'x',xmin-10, tolerance, 1, 0)
    move_attocube(inst, 'y',ymin-10, tolerance, 1, 0)
    move_attocube(inst, 'z',z, tolerance, 1, 0)
    t0=time.time()
    #start the map
    Nx=int((xmax-xmin)/step+1)
    Ny=int((ymax-ymin)/step+1)
    xs=np.linspace(xmin, xmax, Nx)
    ys=np.linspace(ymin, ymax, Ny)
    for i in range(0, len(xs)):
        x=xs[i]
        for j in range(0, len(ys)):
            y=ys[j]
            clear_output(wait=True)
            #move
            move_attocube(inst, 'x',x, tolerance, 0, 0)
            move_attocube(inst, 'y',y, tolerance, 0, 0)
            move_attocube(inst, 'z',z, tolerance, 0, 0)
            time.sleep(waittime)
            
            #read voltage 
            sample = daq.getSample('/%s/demods/0/sample' % device)
            sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
            Rx = sample["x"][0]
            Ry = sample["y"][0]
            R = sample["R"][0]
            
            #read temperature
            T=read_temperature(LS)
            #add to arrays
            xArray = np.append(xArray, Rx)
            yArray = np.append(yArray, Ry)
            RArray = np.append(RArray, R)
            Temp = np.append(Temp, T)
            
            #time
            ts=time.time()-t0
            
            #write to file
            file = open(filename,'a')
            file.write(format(ts, '.15f')+"\t"+format(read_attocube(anc, 'x'), '.15f')+"\t"+format(read_attocube(anc, 'y'), '.15f')+"\t"+format(read_attocube(anc, 'z'), '.15f')+"\t"+format(Rx, '.15f')+'\t'+format(Ry, '.15f')+'\t'+format(R, '.15f')+'\t'+format(T, '.15f')+'\n')
            file.close()
            #plot maps
           # plotMap(xmin, xmax, ymin, ymax, step, xArray, yArray, RArray)
        
        
            print (x, y)
            
def ScanZ (inst,x, y,zmin, zmax, step,tolerance,waittime, filename):
    filename=add_unique_postfix(filename+'.dat')
    t0=time.time()
    
    file = open(filename,'a')
    file.write('Time(s)'+'\t'+'X'+'\t'+'Y'+'\t'+'Z'+'\t'+'Rx'+'\t'+'Ry'+'\t'+'R'+'\t'+'Temp'+'\n')
    file.close()
    timeArray = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    xArray=np.array([])
    yArray=np.array([])
    RArray=np.array([])
    Temp=np.array([])
    
    #take the backlash
    move_attocube(inst, 'z',zmin-10, tolerance, 1, 0)
    move_attocube(inst, 'x',x-10, tolerance, 1, 0)
    move_attocube(inst, 'y',y-10, tolerance, 1, 0)
    
    move_attocube(inst, 'x',x, tolerance, 1, 0)
    move_attocube(inst, 'y',y, tolerance, 1, 0)
    t0=time.time()
    #start the map
    Nz=int((zmax-zmin)/step+1)
   
    zs=np.linspace(zmin, zmax, Nz)
 
    
    for j in range(0, len(zs)):
        z=zs[j]
        clear_output(wait=True)
        #move
        move_attocube(inst, 'x',x, tolerance, 0, 0)
        move_attocube(inst, 'y',y, tolerance, 0, 0)
        move_attocube(inst, 'z',z, tolerance, 0, 0)
        time.sleep(waittime)

        #read voltage 
        sample = daq.getSample('/%s/demods/0/sample' % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        Rx = sample["x"][0]
        Ry = sample["y"][0]
        R = sample["R"][0]

        #read temperature
        T=read_temperature(LS)
        #add to arrays
        xArray = np.append(xArray, Rx)
        yArray = np.append(yArray, Ry)
        RArray = np.append(RArray, R)
        Temp = np.append(Temp, T)

        #time
        ts=time.time()-t0

        #write to file
        file = open(filename,'a')
        file.write(format(ts, '.15f')+"\t"+format(read_attocube(anc, 'x'), '.15f')+"\t"+format(read_attocube(anc, 'y'), '.15f')+"\t"+format(read_attocube(anc, 'z'), '.15f')+"\t"+format(Rx, '.15f')+'\t'+format(Ry, '.15f')+'\t'+format(R, '.15f')+'\t'+format(T, '.15f')+'\n')
        file.close()
        #plot maps
       # plotMap(xmin, xmax, ymin, ymax, step, xArray, yArray, RArray)

def TempRamp(inst, Ti, Tf, rate, filename, plott, printt, cooldown):

    filename=add_unique_postfix(filename+'.dat')
  
    
    file = open(filename,'a')
    #file.write('Time(s)'+'\t'+'T'+'\t'+'X'+'\t'+'Y'+'\t'+'Z'+'\t'+'Rx'+'\t'+'Ry'+'Ry'+'\n')
    file.write('Time(s)'+'\t'+'T'+'\t'+'Rx'+'\t'+'Ry'+'\t'+'R'+'\t'+'theta'+'\t'+"X"+"\t"+"Y"+"\t"+"Z"+'\t'+'I'+'\n')
    file.close()
    
    timeArray = np.array([])
    temperature = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    demod_theta = np.array([])
    demod_I = np.array([])
    posX=np.array([])
    posY=np.array([])
    posZ=np.array([])

    #LS = initialization_lakeshore335(57600)


    T=read_temperature(inst)
    tolerance=0.05;

    #go to initia point 
    
    set_ramp(inst, 1, 1, 0)
    time.sleep(1)
    set_setpoint(inst, 1, Ti)
    time.sleep(1)
    T=read_temperature(inst)
    while abs(T-Ti)>tolerance:
        time.sleep(10)
        T=read_temperature(inst)
        print(T)

    #start ramp 
    set_ramp(inst, 1, 1, rate)
    time.sleep(1)
    set_setpoint(LS, 1, Tf)
    time.sleep(1)

    t0=time.time()
    while  abs(T-Tf)>tolerance:

        T=read_temperature(LS)
        sample = daq.getSample('/%s/demods/0/sample' % device)
        sample2 = daq.getSample('/%s/demods/2/sample' % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        sample["theta"] = np.angle(sample["x"] + 1j * sample["y"])/np.pi*180
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        I = sample2["x"][0]
        theta = sample["theta"][0]
        temperature = np.append(temperature, T)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)
        demod_theta = np.append(demod_theta, theta)
    
        ts=time.time()-t0
        
        posx=read_attocube(anc, 'x')
        posy=read_attocube(anc, 'y')
        posz=read_attocube(anc, 'z')
        #save
        file = open(filename,'a')
        file.write(format(ts, '.15f')+"\t"+format(T, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\t'+format(theta, '.15f')+"\t"+format(posx, '.15f')+"\t"+format(posy, '.15f')+"\t"+format(posz, '.15f')+"\t"+format(I, '.15f')+'\n')
        file.close()

        if plott:
            fig = plt.figure(figsize=(12,12))
            gs = fig.add_gridspec(2, 1)
            ax1 = fig.add_subplot(gs[0, 0])
            #ax2 = fig.add_subplot(gs[1, 0])

            ax1.set_xlabel('temperature (K)')
            ax1.set_ylabel('Demod X(V)')
            #ax2.set_xlabel('temperature (K)')
            #ax2.set_ylabel('Demod theta (deg)')


            ax1.plot(temperature,demod_x,'-o')
            #ax2.plot(temperature,demod_theta,'-o')


            plt.show()
            
            clear_output(wait=True)
        elif printt:
            clear_output(wait=True)
            print(T)
            time.sleep(0.3)
        else:
            time.sleep(0.3)
    if cooldown:
        set_ramp(LS, 1, 1, 0)
        time.sleep(1)
        set_setpoint(LS, 1, Ti)


def coRotHWP(start_pos,end_pos, step_size,offset, filename, plott, printt):
    
    
   

    #Filename & title
    
    filename =add_unique_postfix(filename+'.dat')
    file = open(filename,'a')
    file.write("Time (s)"+"\t"+"Angle (deg)"+"\t"+"Demod x"+'\t'+"Demod y"+'\t'+"R"+"X"+"\t"+"Y"+"\t"+"Z"+"\t"+"T"+'\n')
    file.close()

    #Scan parameter
    scan_range = np.arange(start_pos, end_pos+1, step_size)

    axis1_index = 1        #1:HWP1
    axis2_index = 2        #1: HWP2

    axis2_rot = newport.NewportESP301Axis(controller,axis2_index-1)
    axis2_rot.enable()

    axis1_rot = newport.NewportESP301Axis(controller,axis1_index-1)
    axis1_rot.enable()

    time_constant = 0.3   #unit: second
    go_back = 5          #preventing hysteresis
    real_start = start_pos - go_back


    #Quantities to be measured
    position = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])

    t0=time.time();
    #Scan
    for pos in scan_range:
       # clear_output(wait=True)
        if (pos == start_pos):
            axis2_rot.move(pos-go_back+offset,absolute=True)
            time.sleep(0.5)
            axis1_rot.move((pos-go_back),absolute=True)
        while (axis1_rot.is_motion_done==False): 
            time.sleep(0.5)
            pass
        while (axis2_rot.is_motion_done==False): 
            time.sleep(0.5)
            pass
       


        axis2_rot.move(pos+offset,absolute=True)
        time.sleep(0.5)
        axis1_rot.move(pos,absolute=True)
        
        while (axis1_rot.is_motion_done==False): 
            time.sleep(0.5)
            pass
        while (axis2_rot.is_motion_done==False): 
            time.sleep(0.5)
            pass

        time.sleep(time_constant*4)
        sample = daq.getSample('/%s/demods/0/sample' % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        position = np.append(position, pos)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)

        #Plot
        if plott:
            
            fig = plt.figure(figsize=(5,5))
            gs = fig.add_gridspec(2, 1)
            ax1 = fig.add_subplot(gs[0, 0])
            ax1.set_xlabel('Angle (deg)')
            ax1.set_ylabel('Demod x')
            ax1.plot(position,demod_x,'-o')
            plt.show()
            clear_output(wait=True)
        elif printt:
            clear_output(wait=True)
            print(pos)
            #time.sleep(0.25)
  
        posx=read_attocube(anc, 'x')
        posy=read_attocube(anc, 'y')
        posz=read_attocube(anc, 'z')
        #Save file
        ts=time.time()-t0
        file = open(filename,'a')
        T=read_temperature(LS)
        file.write(format(ts, '.15f')+"\t"+format(pos, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+"\t"+format(posx, '.15f')+"\t"+format(posy, '.15f')+"\t"+format(posz, '.15f')+"\t"+format(T, '.15f')+'\n')
        file.close()
     
        
    angle=-go_back
    axis1_rot.move(angle,absolute=True) 
    angle=-go_back+offset
    axis2_rot.move(angle,absolute=True) 
    while (axis1_rot.is_motion_done==False or axis2_rot.is_motion_done==False):
            pass



#LS=initialization_lakeshore335(57600)

def FreqSweep(freqs, length, fileRoot, plott):
    for f in freqs:
        daq.setDouble('/%s/oscs/0/freq' % device, f)
        print(f)
        fint=int(f)
        fileRoot2=fileRoot+'_%d.dat' %fint
        filename=add_unique_postfix(fileRoot2)
        time.sleep(30)
        t0=time.time()
        file = open(filename,'a')
        file.write('Time(s)'+'\t'+'Demod x'+'\t'+'Demod y'+'\t'+'R'+'\n')
       # file.close()
        timeArray = np.array([])
        demod_x = np.array([])
        demod_y = np.array([])
        demod_r = np.array([])
        demod_theta = np.array([])

        ts=time.time()-t0
        while ts<length:


            sample = daq.getSample('/%s/demods/0/sample' % device)
            sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
            sample["theta"] = np.angle(sample["x"] + 1j * sample["y"])/np.pi*180

            x = sample["x"][0]
            y = sample["y"][0]
            r = sample["R"][0]
            theta = sample["theta"][0]
            ts=time.time()-t0


            timeArray = np.append(timeArray, ts)
            demod_x = np.append(demod_x, x)
            demod_y = np.append(demod_y, y)
            demod_r = np.append(demod_r, r)
            demod_theta = np.append(demod_theta, theta)

            file.write(format(ts, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\n')
            
            if plott:
                clear_output(wait=True)
                fig = plt.figure(figsize=(12,12))
                plt.title('filename')
                gs = fig.add_gridspec(2, 1)
                ax1 = fig.add_subplot(gs[0, 0])
                ax2 = fig.add_subplot(gs[1, 0])
                ax1.set_xlabel('Time (s)')
                ax1.set_ylabel('Demod x')
                ax2.set_xlabel('Time (s)')
                ax2.set_ylabel('theta (deg)')
                ax1.plot(timeArray,demod_x,'-o')
                ax2.plot(timeArray,demod_theta,'-o')

                plt.show()
               # time.sleep(1)
        file.close()
        
def RecordTime(inst, filename,length, plott, printt):

    filename=add_unique_postfix(filename+'.dat')
  
    
    file = open(filename,'a')
    #file.write('Time(s)'+'\t'+'T'+'\t'+'X'+'\t'+'Y'+'\t'+'Z'+'\t'+'Rx'+'\t'+'Ry'+'Ry'+'\n')
    file.write('Time(s)'+'\t'+'T'+'\t'+'Rx'+'\t'+'Ry'+'\t'+'R'+'\t'+'theta'+'\t'+"X"+"\t"+"Y"+"\t"+"Z"+'\n')
    file.close()
    
    timeArray = np.array([])
    temperature = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    demod_r = np.array([])
    demod_theta = np.array([])
    posX=np.array([])
    posY=np.array([])
    posZ=np.array([])

   
    

    t0=time.time()
    ts=time.time()-t0
    while ts<length:
        T=read_temperature(LS)
        sample = daq.getSample('/%s/demods/0/sample' % device)
        sample["R"] = np.abs(sample["x"] + 1j * sample["y"])
        sample["theta"] = np.angle(sample["x"] + 1j * sample["y"])/np.pi*180
        x = sample["x"][0]
        y = sample["y"][0]
        r = sample["R"][0]
        theta = sample["theta"][0]
        ts=time.time()-t0

        temperature = np.append(temperature, T)
        timeArray=np.append(timeArray, ts)
        demod_x = np.append(demod_x, x)
        demod_y = np.append(demod_y, y)
        demod_r = np.append(demod_r, r)
        demod_theta = np.append(demod_theta, theta)



        posx=read_attocube(anc, 'x')
        posy=read_attocube(anc, 'y')
        posz=read_attocube(anc, 'z')
        #save
        file = open(filename,'a')
        file.write(format(ts, '.15f')+"\t"+format(T, '.15f')+"\t"+format(x, '.15f')+'\t'+format(y, '.15f')+'\t'+format(r, '.15f')+'\t'+format(theta, '.15f')+"\t"+format(posx, '.15f')+"\t"+format(posy, '.15f')+"\t"+format(posz, '.15f')+'\n')
        file.close()

        if plott:
            fig = plt.figure(figsize=(12,12))
            gs = fig.add_gridspec(2, 1)
            ax1 = fig.add_subplot(gs[0, 0])
            #ax2 = fig.add_subplot(gs[1, 0])

            ax1.set_xlabel('time(min)')
            ax1.set_ylabel('Demod X(V)')
            #ax2.set_xlabel('temperature (K)')
            #ax2.set_ylabel('Demod theta (deg)')


            ax1.plot(timeArray,demod_x,'-o')
            #ax2.plot(timeArray,demod_theta,'-o')


            plt.show()

            clear_output(wait=True)
        elif printt:
            clear_output(wait=True)
            print(T)
            time.sleep(0.3)
        else:
            time.sleep(0.3)

        
def CorotMapRegZ (xmin, xmax, ymin, ymax,z,step,tolerance,waittime, start_pos,end_pos, step_size,offset,filename, plott, printt):
   
    t0=time.time()
    
    
    timeArray = np.array([])
    demod_x = np.array([])
    demod_y = np.array([])
    xArray=np.array([])
    yArray=np.array([])
    RArray=np.array([])
    
    #take the backlash
    move_attocube(anc, 'x',xmin-10, tolerance, 1, 0)
    move_attocube(anc, 'y',ymin-10, tolerance, 1, 0)
    t0=time.time()
    
    i=0
    j=0
    #start the map
    Nx=int((xmax-xmin)/step+1)
    Ny=int((ymax-ymin)/step+1)
    xs=np.linspace(xmin, xmax, Nx)
    ys=np.linspace(ymin, ymax, Ny)
    for i in range(0, len(xs)):
        x=xs[i]
        for j in range(0, len(ys)):
            y=ys[j]
            clear_output(wait=True)
            #move
            move_attocube(anc, 'x',x, tolerance, 0, 0)
            move_attocube(anc, 'y',y, tolerance, 0, 0)
            move_attocube(anc, 'z',z, tolerance, 0, 0)
            time.sleep(waittime)
            
            filename1=filename+'_'+str(i+1)+'_'+str(j+1)
            
            coRotHWP(start_pos,end_pos, step_size,offset, filename1, plott, printt)

            
            print (x, y)        
