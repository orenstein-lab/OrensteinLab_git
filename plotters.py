"""
Plotting Methods (Plotly Version, Backwards Compatible, Fully Commented)
- Drop-in replacement for matplotlib-based plotting, using Plotly FigureWidget.
- Each variable gets its own subplot (row).
- Tight, readable layout: large margins, rotated x-ticks, font sizes, automargin.
- 1D (line) and 2D (heatmap) plotting.
- Backwards compatible with original API: returns fig, axes (dummy), plot_handles_dict, data containers.
"""

import numpy as np
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from IPython.display import display  # Import display once at the top

# Default figure size: width, height per subplot (in pixels)
SCALE=4/5
DEFAULT_FIGSIZE = (600*SCALE, 300*SCALE)
TEMPLATE = "ggplot2"
FONT = "Arial, sans-serif"
MARKERSIZE = 10
TITLESIZE = 14
TICKSIZE = 12
VSPACE = 0.075/SCALE
EXPONENTFORMATX = 'E'
EXPONENTFORMATY = 'E'
COLORSCALE = 'RdBu'
TITLEDICT =dict(size=TITLESIZE, family=FONT)
TICKDICT =dict(size=TICKSIZE, family=FONT)

def setup_1d_plots_append(vars, xlabel, figsize=DEFAULT_FIGSIZE):
    """
    Set up 1D line plots (one subplot per variable), for data that is appended over time.

    Args:
        vars (list of str): Names of y-axis variables to plot (one subplot per variable).
        xlabel (str): Label for the x-axis (shared across all subplots).
        figsize (tuple): (width, height) in pixels for each subplot (default: (700, 350)).

    Returns:
        fig (go.FigureWidget): The interactive Plotly figure.
        axes (list): Dummy axes (for compatibility with matplotlib API).
        plot_handles_dict (dict): Maps variable name to Plotly trace object.
        xrange (np.ndarray): Array for x-axis data (starts empty).
        vdata1d_dict (dict): Maps variable name to array of y data (starts empty).
    """
    nvars = len(vars)
    vdata1d_dict = {v: np.array([]) for v in vars}
    xrange = np.array([])

    # Create stacked subplots (one per variable), shared x-axis
    fig = make_subplots(
        rows=nvars, cols=1, vertical_spacing=VSPACE
    )
    fig = go.FigureWidget(fig)
    plot_handles_dict = {}

    for i, v in enumerate(vars):
        trace = go.Scatter(x=[], y=[], mode='lines+markers', name=v)
        fig.add_trace(trace, row=i+1, col=1)
        # Set y-axis label for each subplot
        fig.update_yaxes(title_text=v, row=i+1, col=1)
        plot_handles_dict[v] = fig.data[-1]

    # Set x-axis label
    fig.update_xaxes(title_text=xlabel)

    # Tidy layout: increase left margin, rotate x-ticks, set font sizes
    fig.update_layout(
        width=figsize[0],
        height=figsize[1] * nvars,
        showlegend=False,
        template=TEMPLATE,
        margin=dict(l=90, r=30, t=25, b=50),
    )
    fig.update_yaxes(title_font=TITLEDICT, tickfont=TICKDICT, automargin=True, exponentformat=EXPONENTFORMATY)
    fig.update_xaxes(title_font=TITLEDICT, tickfont=TICKDICT, exponentformat=EXPONENTFORMATX)
    fig.update_traces(marker=dict(size=MARKERSIZE))
    display(fig)
    axes = [None] * nvars  # Dummy axes for compatibility
    return fig, axes, plot_handles_dict, xrange, vdata1d_dict

def update_1d_plots_append(fig, axes, vars, plot_handles_dict, xrange, vdata1d_dict, newx, newdata):
    """
    Update 1D plots by appending a new data point to each variable.

    Args:
        fig (go.FigureWidget): The interactive Plotly figure.
        axes (list): Dummy axes (ignored, for compatibility).
        vars (list of str): Names of y-axis variables to update.
        plot_handles_dict (dict): Maps variable name to Plotly trace object.
        xrange (np.ndarray): Current x-axis data.
        vdata1d_dict (dict): Current y data arrays.
        newx (float): New x value to append.
        newdata (dict): Maps variable name to new y value to append.

    Returns:
        xrange (np.ndarray): Updated x-axis data.
        vdata1d_dict (dict): Updated y data arrays.
    """
    xrange = np.append(xrange, newx)
    for v in vars:
        vdata1d_dict[v] = np.append(vdata1d_dict[v], newdata[v])
        trace = plot_handles_dict[v]
        trace.x = xrange
        trace.y = vdata1d_dict[v]
    return xrange, vdata1d_dict

def setup_1d_plots(vars, xlabel, xrange, figsize=DEFAULT_FIGSIZE):
    """
    Set up 1D line plots (one subplot per variable), for data with known x-values.

    Args:
        vars (list of str): Names of y-axis variables to plot (one subplot per variable).
        xlabel (str): Label for the x-axis (shared across all subplots).
        xrange (np.ndarray): Array of x-axis values (predefined).
        figsize (tuple): (width, height) in pixels for each subplot (default: (700, 350)).

    Returns:
        fig (go.FigureWidget): The interactive Plotly figure.
        axes (list): Dummy axes (for compatibility with matplotlib API).
        plot_handles_dict (dict): Maps variable name to Plotly trace object.
        vdata1d_dict (dict): Maps variable name to array of y data (initialized to zeros).
    """
    nvars = len(vars)
    vdata1d_dict = {v: np.zeros(len(xrange)) for v in vars}

    fig = make_subplots(
        rows=nvars, cols=1, vertical_spacing=VSPACE
    )
    fig = go.FigureWidget(fig)
    plot_handles_dict = {}

    for i, v in enumerate(vars):
        trace = go.Scatter(x=xrange, y=vdata1d_dict[v], mode='lines+markers', name=v)
        fig.add_trace(trace, row=i+1, col=1)
        fig.update_yaxes(title_text=v, row=i+1, col=1)
        plot_handles_dict[v] = fig.data[-1]

    fig.update_xaxes(title_text=xlabel)

    fig.update_layout(
        width=figsize[0],
        height=figsize[1] * nvars,
        showlegend=False,
        template=TEMPLATE,
        margin=dict(l=90, r=30, t=25, b=50)
    )
    fig.update_yaxes(title_font=TITLEDICT, tickfont=TICKDICT, automargin=True, exponentformat=EXPONENTFORMATY)
    fig.update_xaxes(title_font=TITLEDICT, tickfont=TICKDICT, exponentformat=EXPONENTFORMATX)
    fig.update_traces(marker=dict(size=MARKERSIZE))
    display(fig)
    axes = [None] * nvars
    return fig, axes, plot_handles_dict, vdata1d_dict

def update_1d_plots(fig, axes, vars, plot_handles_dict, vdata1d_dict, xind, newdata):
    """
    Update 1D plots by modifying the y value at a known x-index for each variable.

    Args:
        fig (go.FigureWidget): The interactive Plotly figure.
        axes (list): Dummy axes (ignored, for compatibility).
        vars (list of str): Names of y-axis variables to update.
        plot_handles_dict (dict): Maps variable name to Plotly trace object.
        vdata1d_dict (dict): Current y data arrays.
        xind (int): Index in the x-array to update.
        newdata (dict): Maps variable name to new y value at index xind.

    Returns:
        vdata1d_dict (dict): Updated y data arrays.
    """
    for v in vars:
        vdata1d_dict[v][xind] = newdata[v]
        trace = plot_handles_dict[v]
        trace.y = vdata1d_dict[v]
    return vdata1d_dict

def setup_2d_plots(vars, xlabel, xrange, ylabel, yrange, figsize=DEFAULT_FIGSIZE):
    """
    Set up 2D heatmap plots (one subplot per variable), for data with known x and y values.

    Args:
        vars (list of str): Names of variables to plot (one subplot per variable).
        xlabel (str): Label for the x-axis (shared across all subplots).
        xrange (np.ndarray): Array of x-axis values.
        ylabel (str): Label for the y-axis (shared across all subplots).
        yrange (np.ndarray): Array of y-axis values.
        figsize (tuple): (width, height) in pixels for each subplot (default: (700, 350)).

    Returns:
        fig (go.FigureWidget): The interactive Plotly figure.
        axes (list): Dummy axes (for compatibility with matplotlib API).
        plot_handles_dict (dict): Maps variable name to Plotly heatmap trace object.
        vdata2d_dict (dict): Maps variable name to 2D np array of data (initialized to zeros).
    """
    nvars = len(vars)
    vdata2d_dict = {v: np.zeros((len(yrange), len(xrange))) for v in vars}

    fig = make_subplots(
        rows=nvars, cols=1, vertical_spacing=VSPACE
    )
    fig = go.FigureWidget(fig)
    plot_handles_dict = {}

    total_subplot_height = 1 - VSPACE*(nvars-1)
    subplot_height = total_subplot_height/nvars
    for i, v in enumerate(vars):
        z = vdata2d_dict[v]
        colorbar_y = subplot_height/2 + i*(subplot_height+VSPACE)
        trace = go.Heatmap(
            z=z,
            x=xrange,
            y=yrange,
            colorbar=dict(title=v, len=0.9/nvars, thickness=18, yanchor="middle", y=colorbar_y, x=1.02, exponentformat=EXPONENTFORMATY),
            colorscale=COLORSCALE,
            zmid=0,
            showscale=True
        )
        fig.add_trace(trace, row=i+1, col=1)
        fig.update_yaxes(title_text=v, row=i+1, col=1)
        plot_handles_dict[v] = fig.data[-1]

    fig.update_xaxes(title_text=xlabel)

    fig.update_layout(
        width=figsize[0],
        height=figsize[1] * nvars,
        showlegend=False,
        template=TEMPLATE,
        margin=dict(l=90, r=30, t=25, b=50)
    )
    fig.update_yaxes(title_font=TITLEDICT, tickfont=TICKDICT, automargin=True, exponentformat=EXPONENTFORMATX)
    fig.update_xaxes(title_font=TITLEDICT, tickfont=TICKDICT, exponentformat=EXPONENTFORMATX)
    display(fig)
    axes = [None] * nvars
    return fig, axes, plot_handles_dict, vdata2d_dict

def update_2d_plots(fig, axes, vars, plot_handles_dict, vdata2d_dict, xind, yind, newdata):
    """
    Update 2D heatmap plots by modifying the value at a given (xind, yind) for each variable.

    Args:
        fig (go.FigureWidget): The interactive Plotly figure.
        axes (list): Dummy axes (ignored, for compatibility).
        vars (list of str): Names of variables to update.
        plot_handles_dict (dict): Maps variable name to Plotly heatmap trace object.
        vdata2d_dict (dict): Current 2D data arrays.
        xind (int): x-index to update.
        yind (int): y-index to update.
        newdata (dict): Maps variable name to new value at [yind, xind].

    Returns:
        vdata2d_dict (dict): Updated 2D data arrays.
    """
    for v in vars:
        vdata2d_dict[v][yind, xind] = newdata[v]
        trace = plot_handles_dict[v]
        trace.z = vdata2d_dict[v]
    return vdata2d_dict
