import numpy as np
from codes.resample_obs_to_finer import resample_obs_to_finer

def rescale(arr_obs_cont_resampled, arr_tem1cont_masked, arr_tem1):
    # Rescaling template
    arr_obs_tem_cont_ratio  = arr_obs_cont_resampled[1] / arr_tem1cont_masked[1]
    arr_obs_tem_cont_ratios = np.concatenate(([arr_tem1cont_masked[0]],
                                              [arr_obs_tem_cont_ratio]), axis=0)
    arr_tem1cont_masked[1]  = arr_tem1cont_masked[1] * arr_obs_tem_cont_ratio
    
    arr_obs_tem_cont_ratios_resampled, arr_tem1_masked = resample_obs_to_finer(arr_obs_tem_cont_ratios, 
                                                                               arr_tem1)
    arr_tem1_masked[1]      = arr_tem1_masked[1]     * arr_obs_tem_cont_ratio
    return arr_obs_cont_resampled, arr_tem1_masked, arr_tem1cont_masked