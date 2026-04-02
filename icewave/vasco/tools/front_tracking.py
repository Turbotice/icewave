import numpy as np

def find_peak_near_x0(signal, x0, window_size=10,subpix=True, subpix_method='gaussian'):
    """
    Find the peak in a 1D signal near a given x0, within a specified window size.

    Parameters:
    signal (numpy array): The 1D signal to search for peaks.
    x0 (int): The x-coordinate around which to search for peaks.
    window_size (int): The size of the window around x0 to search for peaks.

    Returns:
    int: The x-coordinate of the peak found near x0, or None if no peak is found.
    """
    # Define the search range
    start = max(0, x0 - window_size)
    end = min(len(signal), x0 + window_size + 1)

    # Extract the segment of the signal to search for peaks
    segment = signal[start:end]

    # Find the index of the maximum value in the segment
    if len(segment) == 0:
        return None  # No data to search

    peak_index = np.argmax(np.abs(segment))

    # Convert the local index back to the original signal index
    peak_x = start + peak_index

    # subpixel peak detection using quadratic interpolation
    if (peak_index > 0 and peak_index < len(segment) - 1)&(subpix==True):
        x1 = peak_x - 1
        y1 = signal[x1]
        x2 = peak_x
        y2 = signal[x2]
        x3 = peak_x + 1
        y3 = signal[x3]
        if subpix_method=='parabolic':
            A = np.array([[x1**2, x1, 1], [x2**2, x2, 1], [x3**2, x3, 1]])
            B = np.array([y1, y2, y3])
            Coefficients = np.linalg.solve(A, B)
            x_max_subpix = -Coefficients[1] / (2 * Coefficients[0])

        elif subpix_method=='gaussian':
            sub = (np.log(y1)-np.log(y3))/(2*(np.log(y1)+np.log(y3)-2*np.log(y2)))
            x_max_subpix = peak_x + sub
        return x_max_subpix
    elif (peak_index == 0)|(peak_index == len(segment) - 1):
        print("Warning: Peak is at the edge of the search window, subpixel detection may not be accurate.")
        return peak_index
    elif subpix==False:
        return peak_x