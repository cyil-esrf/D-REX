import numpy as np
from scipy.ndimage import median_filter, shift
import h5py
import hdf5plugin
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from gradient_corrections_rx_iron import *

"""
Available functions:

 - tth0_deg = find_exp_tth0(h5file)
 - bkg_mask = background_mask_from_intensity(intensity, percentile=20)
 - fwhm, hmx = fwhm_calculation(x,y) Calculate Full Width at Half Maximum (FWHM)
 - param = fit_strain(bin_centers, distribution)
 - counts, bin_edges, bin_centers = histogram(data, bin_size = 3e-5, limit=5e-4) Calculate histogram of strain distribution
 - strain = strain_tth(data,tth0) Calculate strain from 2theta shift
 - strain = strain_obpitch(data, tth0, obx=250, mainx=5000, zoom=10, correction=True) Calculate strain from 2theta shift with obpitch correction 
"""

#-------------------------------------------------------------
# Find the tth0 based on the experimental set up (position of the detector from the sample)

def find_exp_tth0(h5file):
    with h5py.File(h5file, 'r') as f:
        #read first entry of the h5 file
        h5starter = next(iter(f.keys()))
        mainx = np.abs(np.mean(f[f"{h5starter}/instrument/positioners/mainx"][()]))
        ffz = np.mean(f[f"{h5starter}/instrument/positioners/ffz"][()])
    print(f'mainx = {mainx} and ffz = {ffz}')

    tth0 = np.arctan(ffz/mainx)
    print(f'tth0 values is {tth0} rad ({np.rad2deg(tth0)} deg)')
    return float(np.rad2deg(tth0))

#--------------------------------------------------------------
# Mask background

def background_mask_from_intensity(intensity, percentile=20):
    """
    Background mask from intensity map.
    Pixels below a percentile are considered background.
    """
    vals = intensity[np.isfinite(intensity)]
    thr = np.percentile(vals, percentile)

    bkg_mask = intensity <= thr
    return bkg_mask


#--------------------------------------------------------------
# Calculate FWHM

def lin_interp(x, y, i, half):
    """"
    Linear interpolation to find x value at given y (half max)
    """
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_x(x, y):
    """
    Find x values where y crosses half maximum
    """
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]

def fwhm_calculation(x,y):
    """ 
    Calculate Full Width at Half Maximum (FWHM)
    Inputs:
        x : array-like
            x values (e.g., position, time)
        y : array-like
            y values (e.g., intensity, amplitude)
    Outputs:
        fwhm : float
            Full Width at Half Maximum
        hmx : list
            x values at half maximum
    """
    # find the two crossing points
    hmx = half_max_x(x,y)

    # print the answer
    fwhm = hmx[1] - hmx[0]
    
    return fwhm, hmx

# -------------------------------------------------------------
#Gaussian fit

def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def fit_strain(bin_centers, distribution):
    peak = bin_centers[np.argmax(distribution)]
    p0 = [np.max(distribution), peak, 1e-4]

    params, _ = curve_fit(gaussian, bin_centers, distribution, p0=p0)
    _, x0_fit, _ = params
    return params

# --------------------------------------------------------------
# Histogram of strain distribution
import numpy as np

def histogram(data, bin_size=3e-5, limit=None):
    """
    INPUT
    data : array-like
        Data used to compute histogram.
    bin_size : float
        Width of each histogram bin (must be > 0).
    limit : float or None
        If provided, histogram range is [-limit, limit].
        If None, range is computed from data.

    OUTPUT
    counts : array
        Histogram counts.
    bin_edges : array
        Bin edges.
    bin_centers : array
        Bin centers.
    """

    # Convert to flat numpy array
    data = np.asarray(data).ravel()

    # Keep only finite values (removes NaN and ±inf)
    data = data[np.isfinite(data)]

    if data.size == 0:
        raise ValueError("No finite values in input data.")

    if bin_size <= 0:
        raise ValueError("bin_size must be positive.")

    # ---------------------------------------------------------
    # CASE 1 — automatic limits from data
    # ---------------------------------------------------------
    if limit is None:

        data_min = data.min()
        data_max = data.max()

        if data_max == data_min:
            # all values identical → create one bin
            bins = np.array([data_min, data_min + bin_size])
        else:
            n_bins = int(np.ceil((data_max - data_min) / bin_size))

            # Safety check: avoid accidental huge memory allocation
            if n_bins > 10_000_000:
                raise ValueError(
                    "Too many bins would be created. "
                    "Specify 'limit' or increase bin_size."
                )

            bins = data_min + np.arange(n_bins + 1) * bin_size

        counts, bin_edges = np.histogram(data, bins=bins)

    # ---------------------------------------------------------
    # CASE 2 — user-defined symmetric limit
    # ---------------------------------------------------------
    else:

        if limit <= 0:
            raise ValueError("limit must be positive.")

        bins = np.arange(-limit, limit + bin_size, bin_size)

        mask = (data >= -limit) & (data <= limit)
        counts, bin_edges = np.histogram(data[mask], bins=bins)

    # Compute bin centers
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2

    return counts, bin_edges, bin_centers


# Strain calculation from 2theta shift

def strain_tth(data,tth0):
    diff = (np.radians(data)-np.radians(tth0))
    
    return -0.5 * (diff)/ (np.tan(np.radians(tth0/2)))

def strain_obpitch(data, tth0, obx=250, mainx=5000, zoom=10, correction=True):

    if correction:
        # Apply obpitch correction
        data =  tth_correction(data, tth0, obx, mainx, zoom)
    
    # Flatten and remove NaNs
    valid = data[np.isfinite(data)]


    # Calculate strain from tth0
    strain = strain_tth(data,tth0)
    
    return strain