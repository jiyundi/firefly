import numpy as np
from scipy import interpolate

def resample_obs_to_finer(arr_obs_cont, arr_tem_cont):
    # Resampling -- We want a new obs_flux taken from tem_wave's axis values
    # Resamples 1): make f_obs(x_obs)
    interpo_func_obs_cont = interpolate.interp1d(arr_obs_cont[0], 
                                                 arr_obs_cont[1], 
                                                 kind='linear')
    # Resamples 2): make common wavelength axis w/ same numbers
    mask = (arr_tem_cont[0] >= min(arr_obs_cont[0])) & (arr_tem_cont[0] <= max(arr_obs_cont[0]))
    arr_temobs_common_wave   = arr_tem_cont[0][mask]
    arr_tem_cont_flux_masked = arr_tem_cont[1][mask]
    # Resamples 3): apply f_obs(x) for x in x_tem
    arr_obs_cont_flux_resped = interpo_func_obs_cont(arr_temobs_common_wave)
    # Now, this obs has the same size of tem.
    arr_obs_cont_resampled = np.concatenate(([arr_temobs_common_wave],
                                             [arr_obs_cont_flux_resped]), axis=0)
    arr_tem_cont_masked    = np.concatenate(([arr_temobs_common_wave],
                                             [arr_tem_cont_flux_masked]), axis=0)
    return arr_obs_cont_resampled, arr_tem_cont_masked