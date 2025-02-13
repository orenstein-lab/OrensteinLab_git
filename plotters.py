
'''
Plotting Method
'''
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib
import os
import time
import threading
from tqdm.auto import tqdm
import scipy.optimize as opt
import scipy.interpolate as interp
import inspect
import pickle
from OrensteinLab_git.devices.concurrency_classes import LockedVar, StoppableThread, LockedDict
from OrensteinLab_git.devices import MOTOR_DICT, INSTRUMENT_DICT, META_MOTORS, ACTIVE_MOTORS, ACTIVE_INSTRUMENTS

DEFAULT_FIGSIZE = (6,3)

def setup_1d_plots_append(vars, xlabel, figsize=DEFAULT_FIGSIZE):
    '''
    setup plots where data containers are continuously appended to. Good for situations where one does not know a-priori total number of points to plot

    args:
        - vars:     list of y variable names to plot. creates len(vars) plots
        - xlabel:   name of x axis
    
    returns:
        - fig
        - axes
        - plot_handles_dict:    dictionary containing line handles for each variable in var
        - xrange:               np array for storing x axis data    
        - vdata1d_dict:         dictionary containing data containers for each variable in var
    '''

    # setup data grids
    nvars = len(vars)
    vdata1d_dict = {}
    for v in vars:
        vdata = np.array([])
        vdata1d_dict[v] = vdata
    xrange =  np.array([])

    # setup 1d plot
    fig, axes = plt.subplots(nvars, 1, figsize=(figsize[0],nvars*figsize[1]))
    if type(axes) is matplotlib.axes._axes.Axes:
        axes = [axes]
    for ii, ax in enumerate(axes):
        ax.set_xlabel(xlabel)
        ax.set_ylabel(vars[ii])
        ax.grid(True)
    plot_handles_dict = {}
    for ii, v in enumerate(vars):
        line,  = axes[ii].plot(xrange, vdata1d_dict[v], '-o')
        plot_handles_dict[v] = line
    fig.canvas.draw()
    fig.tight_layout()
    fig.show()

    return fig, axes, plot_handles_dict, xrange, vdata1d_dict

def update_1d_plots_append(fig, axes, vars, plot_handles_dict, xrange, vdata1d_dict, newx, newdata):
    '''
    update plots by appending to data containers

    args:
        - fig
        - axes
        - vars:                 list of y axis variables
        - plot_handles_dict:    dictionary containing line handles for each variable in var
        - xrange:               np array for storing x axis data  
        - vdata1d_dict:         dictionary containing data containers for each variable in var
        - newx:                 value of new x axis position
        - newdata:              dictionary of values to update at index xind for each variable in var

    return:
        - xrange:               modified xrange
        - vdata1d_dict:         modified vdata1d_dict
    '''

    xrange = xrange.append(newx)
    for v in vars:
        line = plot_handles_dict[v]
        vdata1d_dict[v] = vdata1d_dict[v].append(newdata[v])
        vdata = vdata1d_dict[v]
        line.set_data(xrange, vdata)
    for ax in axes:
        ax.relim()
        ax.autoscale()
    fig.canvas.draw()
    fig.canvas.flush_events()

    return xrange, vdata1d_dict

def setup_1d_plots(vars, xlabel, xrange, figsize=DEFAULT_FIGSIZE):
    '''
    setup plots where data container lengths are known

    args:
        - vars:     list of y variable names to plot. creates len(vars) plots
        - xlabel:   name of x axis
        - xrange:   x axis data
    
    returns:
        - fig
        - axes
        - plot_handles_dict:    dictionary containing line handles for each variable in var
        - vdata1d_dict:         dictionary containing data containers for each variable in var
    '''

    # setup data grids
    nvars = len(vars)
    nx = len(xrange)
    vdata1d_dict = {}
    for v in vars:
        vdata = np.zeros(nx)
        vdata1d_dict[v] = vdata

    # setup 1d plot
    fig, axes = plt.subplots(nvars, 1, figsize=(figsize[0],nvars*figsize[1]))
    if type(axes) is matplotlib.axes._axes.Axes:
        axes = [axes]
    for ii, ax in enumerate(axes):
        ax.set_xlabel(xlabel)
        ax.set_ylabel(vars[ii])
        ax.grid(True)
    plot_handles_dict = {}
    for ii, v in enumerate(vars):
        line,  = axes[ii].plot(xrange, vdata1d_dict[v], '-o')
        plot_handles_dict[v] = line
    fig.canvas.draw()
    fig.tight_layout()
    fig.show()

    return fig, axes, plot_handles_dict, vdata1d_dict

def update_1d_plots(fig, axes, vars, plot_handles_dict, vdata1d_dict, xind, newdata):
    '''
    update plots where data container lengths are known

    args:
        - fig
        - axes
        - vars:                 list of y axis variables
        - plot_handles_dict:    dictionary containing line handles for each variable in var
        - vdata1d_dict:         dictionary containing data containers for each variable in var
        - xind:                 index to update for each data container
        - newdata:              dictionary of values to update at index xind for each variable in var

    return:
        - vdata1d_dict:         modified vdata1d_dict
    '''

    for v in vars:
        line = plot_handles_dict[v]
        vdata1d_dict[v][xind] = newdata[v]
        vdata = vdata1d_dict[v]
        line.set_ydata(vdata)
    for ax in axes:
        ax.relim()
        ax.autoscale()
    fig.canvas.draw()
    fig.canvas.flush_events()

    return vdata1d_dict

def setup_2d_plots(vars, xlabel, xrange, ylabel, yrange, figsize=DEFAULT_FIGSIZE):
    '''
    setup 2d plots where data container lengths are known

    args:
        - vars:     list of y variable names to plot. creates len(vars) plots
        - xlabel:   name of x axis
        - xrange:   x axis data
        - ylabel:   name of y axis
        - yrange:   y axis data
    
    returns:
        - fig
        - axes
        - plot_handles_dict:    dictionary containing map handles for each variable in var
        - vdata2d_dict:         dictionary containing data containers for each variable in var
    '''

    # setup data grids
    nvars = len(vars)
    nx = len(xrange)
    ny = len(yrange)
    vdata2d_dict = {}
    for v in vars:
        vdata = np.zeros((ny, nx))
        vdata2d_dict[v] = vdata
    X_coor, Y_coor = np.meshgrid(xrange, yrange)
    extent=[xrange[0], xrange[-1], yrange[0], yrange[-1]]

    # setup 2d plots
    fig, axes = plt.subplots(nvars, 1, figsize=(figsize[0],nvars*figsize[1]))
    if type(axes) is matplotlib.axes._axes.Axes:
        axes = [axes]
    plot_handles_dict = {}
    for ii, v in enumerate(vars):
        vdata = vdata2d_dict[v]
        vax = axes[ii]
        vax.set_xlabel(xlabel)
        vax.set_ylabel(ylabel)
        map = vax.imshow(vdata, cmap="bwr", origin='lower', extent=extent, norm=colors.TwoSlopeNorm(0, vmin=min(vdata.min(), -1e-12), vmax=max(vdata.max(), 1e-12)))
        plot_handles_dict[v] = map
        fig.colorbar(map, ax=vax)
    fig.canvas.draw()
    fig.tight_layout()
    fig.show()

    return fig, axes, plot_handles_dict, vdata2d_dict

def update_2d_plots(fig, axes, vars, plot_handles_dict, vdata2d_dict, xind, yind, newdata):
    '''
    update 2d plots where data container lengths are known

    args:
        - fig
        - axes
        - vars:                 list of variables
        - plot_handles_dict:    dictionary containing map handles for each variable in var
        - vdata1d_dict:         dictionary containing data containers for each variable in var
        - xind:                 x index to update for each data container
        - yind:                 y index to update for each data container
        - newdata:              dictionary of values to update at index [yind,xind] for each variable in var

    return:
        - vdata2d_dict:         modified vdata2d_dict
    '''

    time.sleep(0.01)
    for v in vars:
        map = plot_handles_dict[v]
        vdata2d_dict[v][yind, xind] = newdata[v]
        vdata = vdata2d_dict[v]
        map.set_data(vdata)
        map.set_clim(vmin=min(vdata.min(),0), vmax=max(vdata.max(),0))
    fig.canvas.draw()
    fig.canvas.flush_events()

    return vdata2d_dict