"""
Plotting Methods (Plotly Version, Backwards Compatible)
- Drop-in replacement for matplotlib-based plotting, using Plotly FigureWidget.
- All function signatures and return values match the original API.
"""

import numpy as np
import plotly.graph_objs as go

DEFAULT_FIGSIZE = (600, 300)  # pixels

def setup_1d_plots_append(vars, xlabel, figsize=DEFAULT_FIGSIZE):
    """
    Setup plots where data containers are continuously appended to.
    Returns:
        - fig: Plotly FigureWidget
        - axes: list of None (for compatibility)
        - plot_handles_dict: dict of Plotly traces
        - xrange: np array for x data
        - vdata1d_dict: dict of y data arrays
    """
    nvars = len(vars)
    vdata1d_dict = {v: np.array([]) for v in vars}
    xrange = np.array([])

    fig = go.FigureWidget(layout=go.Layout(
        width=figsize[0], height=figsize[1]*nvars,
        template="plotly_white"
    ))

    plot_handles_dict = {}
    for v in vars:
        trace = go.Scatter(x=[], y=[], mode='lines+markers', name=v)
        fig.add_trace(trace)
        plot_handles_dict[v] = fig.data[-1]

    fig.update_layout(
        xaxis_title=xlabel,
        yaxis_title=", ".join(vars) if nvars == 1 else "",
        showlegend=(nvars > 1)
    )
    fig.show()

    axes = [None]*nvars  # Dummy for backwards compatibility
    return fig, axes, plot_handles_dict, xrange, vdata1d_dict

def update_1d_plots_append(fig, axes, vars, plot_handles_dict, xrange, vdata1d_dict, newx, newdata):
    """
    Update plots by appending to data containers.
    Returns updated xrange and vdata1d_dict.
    """
    xrange = np.append(xrange, newx)
    for v in vars:
        vdata1d_dict[v] = np.append(vdata1d_dict[v], newdata[v])
        trace = plot_handles_dict[v]
        trace.x = xrange
        trace.y = vdata1d_dict[v]
    # axes is a dummy, nothing to update
    return xrange, vdata1d_dict

def setup_1d_plots(vars, xlabel, xrange, figsize=DEFAULT_FIGSIZE):
    """
    Setup plots where data container lengths are known.
    Returns:
        - fig: Plotly FigureWidget
        - axes: list of None (for compatibility)
        - plot_handles_dict: dict of Plotly traces
        - vdata1d_dict: dict of y data arrays
    """
    nvars = len(vars)
    vdata1d_dict = {v: np.zeros(len(xrange)) for v in vars}

    fig = go.FigureWidget(layout=go.Layout(
        width=figsize[0], height=figsize[1]*nvars,
        template="plotly_white"
    ))

    plot_handles_dict = {}
    for v in vars:
        trace = go.Scatter(x=xrange, y=vdata1d_dict[v], mode='lines+markers', name=v)
        fig.add_trace(trace)
        plot_handles_dict[v] = fig.data[-1]

    fig.update_layout(
        xaxis_title=xlabel,
        yaxis_title=", ".join(vars) if nvars == 1 else "",
        showlegend=(nvars > 1)
    )
    fig.show()

    axes = [None]*nvars  # Dummy for backwards compatibility
    return fig, axes, plot_handles_dict, vdata1d_dict

def update_1d_plots(fig, axes, vars, plot_handles_dict, vdata1d_dict, xind, newdata):
    """
    Update plots where data container lengths are known.
    Returns updated vdata1d_dict.
    """
    for v in vars:
        vdata1d_dict[v][xind] = newdata[v]
        trace = plot_handles_dict[v]
        trace.y = vdata1d_dict[v]
    # axes is a dummy, nothing to update
    return vdata1d_dict

def setup_2d_plots(vars, xlabel, xrange, ylabel, yrange, figsize=DEFAULT_FIGSIZE):
    """
    Setup 2D plots where data container lengths are known.
    Returns:
        - fig: Plotly FigureWidget
        - axes: list of None (for compatibility)
        - plot_handles_dict: dict of Heatmap traces
        - vdata2d_dict: dict of 2D np arrays
    """
    nvars = len(vars)
    vdata2d_dict = {v: np.zeros((len(yrange), len(xrange))) for v in vars}

    fig = go.FigureWidget(layout=go.Layout(
        width=figsize[0], height=figsize[1]*nvars,
        template="plotly_white"
    ))

    plot_handles_dict = {}
    for v in vars:
        z = vdata2d_dict[v]
        trace = go.Heatmap(
            z=z, x=xrange, y=yrange,
            colorbar=dict(title=v),
            colorscale='RdBu', zmid=0
        )
        fig.add_trace(trace)
        plot_handles_dict[v] = fig.data[-1]

    fig.update_layout(
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        showlegend=False
    )
    fig.show()

    axes = [None]*nvars  # Dummy for backwards compatibility
    return fig, axes, plot_handles_dict, vdata2d_dict

def update_2d_plots(fig, axes, vars, plot_handles_dict, vdata2d_dict, xind, yind, newdata):
    """
    Update 2D plots where data container lengths are known.
    Returns updated vdata2d_dict.
    """
    for v in vars:
        vdata2d_dict[v][yind, xind] = newdata[v]
        trace = plot_handles_dict[v]
        trace.z = vdata2d_dict[v]
    # axes is a dummy, nothing to update
    return vdata2d_dict
