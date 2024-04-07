import numpy as np
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from specutils import Spectrum1D
from specutils.analysis import correlation
from scipy import signal

from codes.find_sus_O_II_line    import find_sus_O_II_line
from codes.resample_obs_to_finer import resample_obs_to_finer
from codes.rescale               import rescale

def cross_corr(arr_obs, arr_tem1,
               z_guess1,
               dic_emi_lines):
    """
    FLUX-CONTINUUM NORMALIZATION AND SUBSTRACTION
    
    # 1) Estimate continuum by ignoring all obscured lines
         "kernel_size=int((len(arr_obs_flux)/5)//2*2 + 1"
         "//2*2 + 1": kernal_size must be odd. Do not change.
                 "5": self-def. Smaller -> fewer obscured lines
               
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
    
    # Observed spec (for all spec1d figures)
    arr_obs_cont_flux = signal.medfilt(arr_obs[1], 
                                       kernel_size=int((len(arr_obs[1])/5)//2*2 + 1)
                                       ) # find continuum
    arr_obs_cont = np.concatenate(([arr_obs[0]],
                                   [arr_obs_cont_flux],
                                   [arr_obs[2]]), axis=0)
    arr_obs_rela_flux = arr_obs[1] / arr_obs_cont[1] - 1 # normalize and subtract
    arr_obs_rela_erro = arr_obs[2] / arr_obs_cont[1] - 1
    arr_obs_rela = np.concatenate(([arr_obs[0]],
                                   [arr_obs_rela_flux],
                                   [arr_obs_rela_erro]), axis=0)
    del arr_obs_rela_flux, arr_obs_rela_erro, arr_obs_cont_flux
    
    # Template-1 spec (to be redshifted, for Figure 1, 2, 3, 4)
    arr_tem1cont_flux = signal.medfilt(arr_tem1[1], 
                                       kernel_size=int((len(arr_obs[1])/5)//2*2 + 1)
                                       )
    arr_tem1cont = np.concatenate(([arr_tem1[0]],
                                   [arr_tem1cont_flux],
                                   [arr_tem1[2]]), axis=0)
    arr_tem1rela_flux = arr_tem1[1] / arr_tem1cont[1] - 1
    arr_tem1rela_erro = arr_tem1[2] / arr_tem1cont[1] - 1
    arr_tem1rela = np.concatenate(([arr_tem1[0]],
                                   [arr_tem1rela_flux],
                                   [arr_tem1rela_erro]), axis=0)
    del arr_tem1cont_flux, arr_tem1rela_flux, arr_tem1rela_erro
    
    # Find z
    # fit
    obs_rela_spec = Spectrum1D(spectral_axis = arr_obs_rela[0] * u.AA, 
                               flux          = arr_obs_rela[1] * u.dimensionless_unscaled, 
                               uncertainty   = StdDevUncertainty(
                                               arr_obs_rela[2] * u.dimensionless_unscaled))
    tem1rela_spec = Spectrum1D(spectral_axis = arr_tem1rela[0] * u.AA, 
                               flux          = arr_tem1rela[1] * u.dimensionless_unscaled, 
                               uncertainty   = StdDevUncertainty(
                                               arr_tem1rela[2] * u.dimensionless_unscaled))
    corr1, lag1   = correlation.template_correlate(observed_spectrum = obs_rela_spec, 
                                                   template_spectrum = tem1rela_spec,
                                                   lag_units=u.dimensionless_unscaled)
    corr1 = corr1[lag1 > 0.]
    lag1  =  lag1[lag1 > 0.]
    corr1_normalized = corr1/sum(corr1)
    
    # find peak z
    z_peak1 = lag1[np.argmax(corr1_normalized)].value
    
    # no action if z is otherwise specified
    if z_guess1 == None:
        z_guess1 = z_peak1

    # redshift tem spec by (1+z)
    arr_tem1C     = arr_tem1.copy()    # backups
    arr_tem1contC = arr_tem1cont.copy()
    arr_tem1relaC = arr_tem1rela.copy()
    arr_tem1[0]     *= (1+z_guess1)
    arr_tem1cont[0] *= (1+z_guess1)
    arr_tem1rela[0] *= (1+z_guess1)
    
    # Template-1 spec (redshifted + rescaled, for Figure 1, 2)
    arr_obs_cont_resampled, arr_tem1cont_masked = resample_obs_to_finer(arr_obs_cont, 
                                                                        arr_tem1cont)
    # rescale
    arr_obs_cont_resampled, arr_tem1_masked, arr_tem1cont_masked = rescale(arr_obs_cont_resampled,
                                                                           arr_tem1cont_masked, 
                                                                           arr_tem1)

    z_guess1, (arr_tem1_masked_newz,     # Figure 3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
               arr_tem1cont_masked_newz, # Figure 3, 4/4
               arr_tem1relaC             # Figure 4, 2/2
               ) = find_sus_O_II_line(z_guess1, z_peak1,
                                      arr_obs,
                                      arr_obs_cont,
                                      arr_obs_rela, 
                                      arr_tem1C,
                                      arr_tem1contC,
                                      arr_tem1relaC,
                                      dic_emi_lines)

    # z_guess2, result_dic4 = find_sus_O_II_line(z_guess2, z_peak2,
    #                                            arr_obs_wave, arr_obs_flux, arr_obs_erro,
    #                                            arr_obs_rela_flux, 
    #                                            arr_tem2wave, arr_tem2flux, arr_tem2erro,
    #                                            dic_emi_lines)
    # arr_tem2wave_z2          = result_dic3['arr_tem_wave_z']
    # arr_tem2rela_flux_z2     = result_dic3['arr_tem_rela_flux_z']
    # arr_tem2wave_new_z2      = result_dic3['arr_tem_wave_new_z']
    # arr_tem2flux_new_z2      = result_dic3['arr_tem_flux_new_z']
    # arr_tem2cont_flux_new_z2 = result_dic3['arr_tem_cont_flux_new_z']

    
    # result_dic = {
    #     # 'z_guess': z_guess,
    #     # 'z_peak': z_peak,
    #     # 'lag': lag,
    #     # 'corr_normalized': corr_normalized,
    #     # 'arr_obs_rela_flux': arr_obs_rela_flux,
    #     # 'arr_obs_cont_wave': arr_obs_cont_wave,
    #     # 'arr_obs_cont_flux': arr_obs_cont_flux,
    #     # 'arr_tem_rela_flux': arr_tem_rela_flux,
    #     # 'arr_tem_wave_new': arr_tem_wave_new,
    #     # 'arr_tem_flux_new': arr_tem_flux_new,
    #     # 'arr_tem_cont_flux_new': arr_tem_cont_flux_new
    # }
    return (z_peak1, z_guess1,
    lag1, corr1_normalized, 
    arr_obs_cont,        # Figure 1+3: 2/4 (1/4=arr_obs)
    arr_tem1_masked,     # Figure 1: 3/4 (Can't use arr_tem1 b/c un-normalized!)
    arr_tem1cont_masked, # Figure 1, 4/4
    arr_obs_rela,        # Figure 2+4, 1/2
    arr_tem1rela,        # Figure 2, 2/2
    arr_tem1_masked_newz,     # Figure 3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
    arr_tem1cont_masked_newz, # Figure 3, 4/4
    arr_tem1relaC             # Figure 4, 2/2
    )