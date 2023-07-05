import instruments.newport as newport
import os.path
import time

#ESP301 initialization
f_conf = open(os.path.dirname(__file__)+ r'\..\..\Configuration.txt', "r")
conf_info = f_conf.read()
conf_info_split = conf_info.split('\n')
port_id = conf_info_split[1].split('\t')[1]
f_conf.close()
controller = newport.NewportESP301.open_serial(port=port_id, baud=921600)
time.sleep(0.1)
axes = [1, 2, 3]
for axis in axes:
    axis_obj = newport.NewportESP301Axis(controller, axis - 1)
    axis_obj.enable()
    time.sleep(0.1)

def read(axis):
    axis_obj = newport.NewportESP301Axis(controller, axis - 1)
    return axis_obj.position

def move_absolute(axis, value):
    axis_obj = newport.NewportESP301Axis(controller, axis - 1)
    axis_obj.move(value, absolute=True)

def move_relative(axis, value):
    axis_obj = newport.NewportESP301Axis(controller, axis - 1)
    axis_obj.move(value, absolute=False)

def stab_absolute(axis, value):
    axis_obj = newport.NewportESP301Axis(controller, axis - 1)
    axis_obj.move(value, absolute=True)
    while axis_obj.is_motion_done == False:
        time.sleep(0.1)

def stab_absolute_duo(axis1, value1, axis2, value2):
    axis_obj1 = newport.NewportESP301Axis(controller, axis1 - 1)
    axis_obj1.move(value1, absolute=True)
    time.sleep(0.1)
    axis_obj2 = newport.NewportESP301Axis(controller, axis2 - 1)
    axis_obj2.move(value2, absolute=True)
    time.sleep(0.1)
    while axis_obj1.is_motion_done == False or axis_obj2.is_motion_done == False:
        time.sleep(0.1)

def stab_relative(axis, value):
    axis_obj = newport.NewportESP301Axis(controller, axis - 1)
    axis_obj.move(value, absolute=False)
    while axis_obj.is_motion_done == False:
        time.sleep(0.1)

def isMoving(axis):
    axis_obj = newport.NewportESP301Axis(controller, axis - 1)
    try:
        motion_boolean = axis_obj.is_motion_done
    except ValueError:
        return False
    return motion_boolean