# -*- coding: utf-8 -*-
"""
Created on Wed Apr  9 16:39:25 2025

@author: sebas
"""
import numpy as np
from matplotlib.path import Path
import matplotlib.pyplot as plt


def interactive_scatter_filter(x, y,properties = None):
    """
    Launches an interactive matplotlib window to manually filter points from a scatter plot
    by drawing polygons.

    Parameters:
        x (array-like, len N): x coordinates of the scatter points.
        y (array-like, len N): y coordinates of the scatter points.
        properties (array-like, N x d) or dictionnary containing 1D-field of length N : properties of the plotted points
            default is None

    Returns:
        x_filtered (np.ndarray): Filtered x coordinates.
        y_filtered (np.ndarray): Filtered y coordinates.
        filtered_properties (same object as input): filtered properties of plotted points
    """
    points = np.column_stack((x, y))
    filtered_points = points.copy()
    filtered_properties = properties.copy()
    polygon_vertices = []

    fig, ax = plt.subplots()
    ax.scatter(filtered_points[:, 0], filtered_points[:, 1])
    plt.title("Left-click: add point | Right/Middle-click: close polygon | 'q': quit")

    # Use a flag to track when the user is done
    done = {'value': False}
    
    def redraw():
        ax.clear()
        ax.scatter(filtered_points[:, 0], filtered_points[:, 1])
        ax.set_title("Left-click: add point | Right/Middle-click: close polygon | 'q': quit")
        fig.canvas.draw()

    def on_click(event):
        nonlocal filtered_points, polygon_vertices, filtered_properties
        
        # if event.key == 'q':
        #     print("Exiting polygon selection.")
        #     fig.canvas.mpl_disconnect(cid_click)
        #     fig.canvas.mpl_disconnect(cid_key)
        #     return
    
        if event.button == 1 and event.inaxes:  # Left click
            polygon_vertices.append((event.xdata, event.ydata))
            ax.plot(event.xdata, event.ydata, 'ro')
            fig.canvas.draw()

        elif event.button in [2, 3] and len(polygon_vertices) >= 3:  # Close polygon
            polygon = Path(polygon_vertices)
            inside = polygon.contains_points(filtered_points)
            filtered_points = filtered_points[~inside]
            
            if filtered_properties is not None :
                if type(filtered_properties) == dict:
                    for key in filtered_properties.keys():
                        filtered_properties[key] = filtered_properties[key][~inside]
                elif type(filtered_properties) == np.ndarray :
                    for i in range(filtered_properties.shape[-1]):
                        filtered_properties[:,i] = filtered_properties[~inside,i]
                else :
                    print('Properties are not updated, wrong format !')
                    
            polygon_vertices = []
            redraw()

    def on_key(event):
        if event.key == 'q':
            print("Exiting polygon selection.")
            done['value'] = True
            fig.canvas.mpl_disconnect(cid_click)
            fig.canvas.mpl_disconnect(cid_key)
            plt.close(fig)

    cid_click = fig.canvas.mpl_connect('button_press_event', on_click)
    cid_key = fig.canvas.mpl_connect('key_press_event', on_key)

    while not done['value']:
        plt.pause(0.1)  # Allows event loop to run

    return filtered_points[:, 0], filtered_points[:, 1], filtered_properties