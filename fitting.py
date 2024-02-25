import numpy as np
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from specutils import Spectrum1D
from specutils.analysis import correlation
from scipy import interpolate
from scipy import signal

def cross_corr(arr_obs_wave, arr_obs_flux, arr_obs_erro,
               arr_tem_wave, arr_tem_flux, arr_tem_erro,
               z_guess):
    """
    FLUX-CONTINUUM NORMALIZATION AND SUBSTRACTION
    
    # 1) Estimate continuum by ignoring all obscured lines
         "//2*2 + 1": kernal_size must be odd
                 "9": self-def. Smaller -> fewer obscured lines
               
    # 2) Calculate relative (template-continuum / obs-flux-continuum)
    
    # 3) Since the correlation function gives a high peak when observed
         and template flux values are both positive or negative. 
         The relative flux in step 2) is always positive. Now subtract 
         the relative by 1, so this makes correlation well-behaved.
    
    # 4) Fitting -
         Assemble [wavelength, spec, error] into the spec objects
         spectral_axis: array of float, quantity with unit u.AA
         flux:          array of float, quantity with unit u.Jy, 
                        u.Unit("erg / (Angstrom s cm2)"), or defined
         uncertainty:   array of float, quantity with unit u.Jy, 
                        u.Unit("erg / (Angstrom s cm2)"), or defined
         This step gives you a redshift probability distribution. 
         Only redshifts > 0 are realistic.
    
    # 5) Resampling and Rescaling by continuum
    """
    arr_obs_cont_flux = signal.medfilt(arr_obs_flux, 
                                       kernel_size=int((len(arr_obs_flux)/9)//2*2 + 1)
                                       )
    arr_tem_cont_flux = signal.medfilt(arr_tem_flux, 
                                       kernel_size=int((len(arr_obs_flux)/9)//2*2 + 1)
                                       )
    
    arr_obs_rela_flux = arr_obs_flux / arr_obs_cont_flux - 1
    arr_obs_rela_erro = arr_obs_erro / arr_obs_cont_flux - 1
    arr_tem_rela_flux = arr_tem_flux / arr_tem_cont_flux - 1
    arr_tem_rela_erro = arr_tem_erro / arr_tem_cont_flux - 1
    
    obs_rela_spec = Spectrum1D(spectral_axis = arr_obs_wave     , 
                               flux          = arr_obs_rela_flux, 
                               uncertainty   = StdDevUncertainty(
                                               arr_obs_rela_erro))
    tem_rela_spec = Spectrum1D(spectral_axis = arr_tem_wave     , 
                               flux          = arr_tem_rela_flux, 
                               uncertainty   = StdDevUncertainty(
                                               arr_tem_rela_erro))
    corr, lag = correlation.template_correlate(observed_spectrum = obs_rela_spec, 
                                               template_spectrum = tem_rela_spec,
                                               lag_units=u.dimensionless_unscaled)
    corr = corr[lag > 0.]
    lag  =  lag[lag > 0.]
    corr_normalized = corr/sum(corr)
    z_peak = lag[np.argmax(corr_normalized)].value
    if z_guess == None:
        z_guess = z_peak

    arr_obs_cont_wave = arr_obs_wave     
    arr_tem_wave      *= (1+z_guess)
    
    # Resampling -- We want a new obs_flux taken from tem_wave's axis values
    # Resamples 1): make f_obs(x_obs)
    interpo_func_obs_cont = interpolate.interp1d(arr_obs_wave     , 
                                                 arr_obs_cont_flux, 
                                                 kind='linear')
    # Resamples 2): make common wavelength axis w/ same numbers
    mask = (arr_tem_wave >= min(arr_obs_wave)) & (arr_tem_wave <= max(arr_obs_wave))
    arr_temobs_common_wave   =      arr_tem_wave[mask]
    arr_tem_cont_flux_new    = arr_tem_cont_flux[mask]
    # Resamples 3): apply f_obs(x) for x in x_tem
    arr_obs_cont_flux_new    = interpo_func_obs_cont(arr_temobs_common_wave)
    # Rescaling template
    arr_obs_tem_cont_ratio  = arr_obs_cont_flux_new / arr_tem_cont_flux_new
    arr_tem_cont_flux_new   = arr_tem_cont_flux[mask] * arr_obs_tem_cont_ratio
    arr_tem_flux_new        = arr_tem_flux[mask] * arr_obs_tem_cont_ratio
    arr_tem_wave_new        = arr_tem_wave[mask]
    
    result_dic = {
        'z_guess': z_guess,
        'z_peak': z_peak,
        'lag': lag,
        'corr_normalized': corr_normalized,
        'arr_obs_rela_flux': arr_obs_rela_flux,
        'arr_obs_cont_wave': arr_obs_cont_wave,
        'arr_obs_cont_flux': arr_obs_cont_flux,
        'arr_tem_rela_flux': arr_tem_rela_flux,
        'arr_tem_wave_new': arr_tem_wave_new,
        'arr_tem_flux_new': arr_tem_flux_new,
        'arr_tem_cont_flux_new': arr_tem_cont_flux_new
    }
    return result_dic