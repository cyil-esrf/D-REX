import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.signal import detrend
from scipy.io import savemat

'''
Available functions:
- chi_corr = chi_com_correction(chi_cm, obpitch, obx = 250 , mainx = 5000, zoom = 10, doRot90=False)
- twotheta_corr_deg =  tth_correction(twotheta_meas_deg,obpitch, obx = 250, mainx = 5000, zoom = 10) 
'''

#-------------------------------------------------------------
#                     chi correction
#-------------------------------------------------------------


def chi_com_correction(chi_cm, obpitch, obx = 250 , mainx = 5000, zoom = 10, doRot90=False):

    # --- unit conversion ---
    obx = obx * 1e-3     # mm → m
    mainx = mainx * 1e-3
    
    molt = zoom/10
    
    d1=obx/np.cos(np.deg2rad(obpitch))
    d2=(np.abs(mainx)/np.cos(np.deg2rad(obpitch))) - d1
    pix=40e-9 / molt
    N=88
    R=50e-6
    rho=1.845
    delta_over_rho=6.3930e-7
    mu_over_rho=0.2558
    twoThetaB_deg=obpitch
    E_eV = 17000

    """
    Apply chi-COM optics correction (Poulsen et al., Eq. 27)
    and zero-mean normalization.

    Parameters
    ----------
    chi_cm: chi data np.array
    doRot90 : bool
        Rotate chi map by 90 deg (to match workflow)
    d1, d2 : float
        Distances [m]
    pix : float
        Pixel size [m]
    N, R : float
        CRL parameters
    rho, delta_over_rho, mu_over_rho : float
        Material parameters
    twoThetaB_deg : float
        Bragg angle (2θ_B) in degrees

    Returns chi_corr_zero : corrected chi-COM
    """

    if doRot90:
        chi_cm = np.rot90(chi_cm)

    nr, nc = chi_cm.shape

    # -------------------- GEOMETRY -----------------------------
    M = d2 / d1
    theta = np.deg2rad(twoThetaB_deg / 2)

    # -------------------- MATERIAL -----------------------------
    delta = delta_over_rho * rho
    mu = mu_over_rho * rho * 100.0   # cm^-1 → m^-1

    # -------------------- REAL-SPACE y_s ------------------------
    # CORRECT CENTERING (Python, 0-based indexing)
    CC = np.arange(nc)[None, :]
    y_s = (CC - nc / 2) * pix
    # y_s = np.repeat(y_s, nr, axis=0)

    # -------------------- sigma_a (Eq. 9) ----------------------
    sigma_a = delta * (M / (M + 1)) * np.sqrt((2 * N) / (mu * R))
    NA_FWHM = 2.35 * sigma_a

    # -------------------- xi (large d1/f', weak vignetting) -----
    xi = 1.0 / d1

    # -------------------- APPLY EQ. (27) -----------------------
    eta_shift_rad = (xi * y_s) / (2 * np.sin(theta))
    eta_shift_deg = np.rad2deg(eta_shift_rad)

    chi_corr = chi_cm - eta_shift_deg

    # -------------------- ZERO-MEAN ----------------------------
    # chi_corr_zero = chi_corr - np.nanmean(chi_corr)


    # -------------------- RETURN -------------------------------
    # return chi_corr_zero
    return chi_corr



#-------------------------------------------------------------
#                   obpitch correction
#-------------------------------------------------------------


def tth_correction(twotheta_meas_deg,obpitch, obx = 250, mainx = 5000, zoom = 10):

    # --- unit conversion ---
    obx = obx * 1e-3     # mm → m
    mainx = mainx * 1e-3
    
    molt = zoom/10
    
    d1=obx/np.cos(np.deg2rad(obpitch))
    d2=(np.abs(mainx)/np.cos(np.deg2rad(obpitch))) - d1
    pix=40e-9 / molt
    N=88
    R=50e-6
    rho=1.845
    delta_over_rho=6.3930e-7
    mu_over_rho=0.2558
    twoThetaB_deg=obpitch
    
    """
    Apply Eq.(28) 2θ COM correction and compute strain.

    Parameters
    ----------
    twotheta_meas_deg : obpitch data np.array
    d1, d2 : float
        Distances defining magnification
    pix : float
        Effective pixel size [m]
    twoThetaB_deg : float
        Bragg angle 2θ_B [deg]
    doRot90 : bool
        Rotate input map by 90° (MATLAB compatibility)
    make_plots : bool
        Show diagnostic plots
    save_npz : bool
        Save results to NPZ
    outname : str
        Output filename

    Returns
    -------
    results : dict
        Dictionary with corrected 2θ and strain fields
    """

    # -------------------- GEOMETRY -------------------------
    M = d2 / d1
    xi = 1 / d1                      # 1/m


    nr, nc = twotheta_meas_deg.shape

    # -------------------- z_s (CORRECT) --------------------
    # MATLAB: (RR - (nr+1)/2) * pix
    z_s = -(np.arange(nr)[:, None] - nr / 2) * pix

    # -------------------- Eq. (28) -------------------------
    twotheta_meas_rad = np.deg2rad(twotheta_meas_deg)
    twotheta_shift_rad = -xi * z_s
    twotheta_corr_rad = twotheta_meas_rad - twotheta_shift_rad
    twotheta_corr_deg = np.rad2deg(twotheta_corr_rad)

    # # -------------------- DETREND --------------------------
    # x = np.arange(nr)
    # A = np.vstack([x, np.ones(nr)]).T

    # twotheta_detr_deg = np.empty_like(twotheta_meas_deg)
    # for j in range(nc):
    #     m, c = np.linalg.lstsq(A, twotheta_meas_deg[:, j], rcond=None)[0]
    #     twotheta_detr_deg[:, j] = twotheta_meas_deg[:, j] - (m * x + c)

    
    
    return twotheta_corr_deg
