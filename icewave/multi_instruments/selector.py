import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import RectangleSelector

import icewave.multi_instruments.situation_map as maps

sBox = None

# Function to handle the selection
def select1(eclick, erelease):
    global sBox1
    # Get the coordinates of the selected region
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print(f'Selected region: x1={x1}, y1={y1}, x2={x2}, y2={y2}')
    sBox1 = [x1,x2,y1,y2]

def select2(eclick, erelease):
    global sBox2
    # Get the coordinates of the selected region
    x1, y1 = eclick.xdata, eclick.ydata
    x2, y2 = erelease.xdata, erelease.ydata
    print(f'Selected region: x1={x1}, y1={y1}, x2={x2}, y2={y2}')
    sBox2 = [x1,x2,y1,y2]

def selector(fig,ax):
    # Create a RectangleSelector
    rectangle_selector = RectangleSelector(ax, on_select,
                                       useblit=True,
                                       button=[1],  # Left mouse button
                                       minspanx=5, minspany=5,
                                       spancoords='pixels',
                                       interactive=True)
    return sBox

def toggle_selector(event):
    print('Key pressed.')
    if event.key == 't':
        for selector in selectors:
            name = type(selector).__name__
            if selector.active:
                print(f'{name} deactivated.')
                selector.set_active(False)
            else:
                print(f'{name} activated.')
                selector.set_active(True)

def main():
    date = '0226'
    disk = 'Fabien_2024/Share_hublot'
    records = maps.get_record(date,year='2024',disk=disk)
    #maps.display_map
    fig,axmap,axtime = maps.synthesis(records,project=False)#,ax=None,BBox=None,project=True)
    # Show the plot
    axs = [axmap,axtime]

    selectors=[]
    for ax, selector_class,on_select in zip(axs, [RectangleSelector, RectangleSelector],[select1,select2]):
        ax.set_title(f"Click and drag to draw a {selector_class.__name__}.")

        selectors.append(selector_class(
            ax,on_select,
            useblit=True,
            button=[1, 3],  ## disable middle button
            minspanx=5, minspany=5,
            spancoords='pixels',
            interactive=True))
        fig.canvas.mpl_connect('key_press_event', toggle_selector)
    #sBox = selector(fig,axmap)
    
#    rectangle_selector.set_active(True)
    #tBox = selector(fig,axtime)

    #print(tBox)
    #print(sBox,tBox)
#    fig.canvas.mpl_connect('key_press_event', toggle_selector)
    plt.show()
    print(sBox1,sBox2)

    space_box = sBox1
    time_box = sBox2

    tmin = time_box[0]
    tmax = time_box[1]

    print(space_box)
    print(tmin,tmax)

    return space_box,tmin,tmax
    

if __name__ == '__main__':
    main()
