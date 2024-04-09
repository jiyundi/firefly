import numpy as np
from scipy.signal import find_peaks

from codes.resample_obs_to_finer import resample_obs_to_finer
from codes.rescale               import rescale

def find_sus_O_II_line(z_guess, z_peak,
                       arr_obs,
                       arr_obs_cont,
                       arr_obs_rela, 
                       arr_tem1C,
                       arr_tem1contC,
                       arr_tem1relaC,
                       dic_emi_lines):
    z_2ndbest = None
    if z_guess == z_peak:
        # Find neighberhood lines and spacings
        i_peaks, prop_dic = find_peaks(arr_obs_rela[1], 
                                       height=float(0.2*np.max(arr_obs_rela[1])))
        w_peaks      = arr_obs[0][i_peaks]
        if len(w_peaks) != 0:
            w_diff_peaks = np.append([w_peaks[0] - arr_obs[0][0]], np.diff(w_peaks))
            w_diff_peaks = np.append(w_diff_peaks, [arr_obs[0][-1] - w_peaks[-1]])
            # i_PEAK0            = i_peaks[arr_obs_rela_flux[i_peaks].argmax()]
            # i_PEAK0_in_w_peaks = np.where(i_peaks == i_PEAK0)[0]
            # In nearby 1200A region, if only one emission line, then it's likely a [O II]
            for i in range(len(w_peaks)):
                i_peak = i_peaks[i]
                w_diff_value_left  = w_diff_peaks[i]
                w_diff_value_right = w_diff_peaks[i+1]
                if w_diff_value_left > 600 and w_diff_value_right > 600:
                    z_2ndbest = arr_obs[0][i_peak] / dic_emi_lines['[OII]'][0] - 1
                    
                    # copied from cross_corr in fitting.py ----------------------------------------------------------------------------
                    arr_tem1C[0]     = arr_tem1C[0]     * (1+z_2ndbest)
                    arr_tem1contC[0] = arr_tem1contC[0] * (1+z_2ndbest)
                    arr_tem1relaC[0] = arr_tem1relaC[0] * (1+z_2ndbest)
                    arr_obs_cont_resampled_newz, arr_tem1cont_masked_newz = resample_obs_to_finer(arr_obs_cont, 
                                                                                                  arr_tem1contC)
                    # rescale
                    arr_obs_cont_resampled_newz, arr_tem1_masked_newz, arr_tem1cont_masked_newz = rescale(arr_obs_cont_resampled_newz,
                                                                                                          arr_tem1cont_masked_newz, 
                                                                                                          arr_tem1C)
                    results = (arr_tem1_masked_newz,     # Figure 3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
                               arr_tem1cont_masked_newz, # Figure 3, 4/4
                               arr_tem1relaC             # Figure 4, 2/2
                               )
                    break
        if z_2ndbest==None:
            z_guess = z_peak
            results = ([[], [], []], 
                       [[], [], []], 
                       [[], [], []])
        else:
            z_guess = z_2ndbest
    return z_guess, results