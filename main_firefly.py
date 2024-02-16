import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.nddata import StdDevUncertainty
from specutils import Spectrum1D
from specutils.analysis import correlation
import sys
import os

import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import signal
from scipy.optimize import curve_fit

from binning_1d     import binning_1d
from load_templates import read_templates
from load_lines     import read_lines
from plot_lines     import plot_lines
global dic_templates, dic_emi_lines, dic_abs_lines
dic_templates                = read_templates()
dic_emi_lines, dic_abs_lines = read_lines()
print('Template and lines loading complete.')
binfactor = 24
z_guess   = None

# USER INPUTS...
inputlist = sys.argv

if len(inputlist) > 1:
    spec1dfolderpath = inputlist[1]
    spec1dobjnames   = inputlist[2]
    templatename     = inputlist[3]
    
    # check illegal inputs
    if spec1dfolderpath[-1] != '/':
        spec1dfolderpath += '/'
        print(    'Spec1d Folder:   ', spec1dfolderpath, 'is specified.')
        
    if spec1dobjnames == '-all':
        input_spec_list = []
        short_list      = os.listdir(spec1dfolderpath)
        for f in short_list: input_spec_list += [spec1dfolderpath + f]
        print(    'File selected:   ', len(input_spec_list), 'files in this spec1d folder.')
    else:
        input_spec_list = [spec1dfolderpath + spec1dobjnames]
        if len(inputlist) <= 4:
            print('File selected:   ', input_spec_list, '(only 1 file).',
                  'You had no redshift guess specified for this object.')
        else:
            z_guess = float(inputlist[4])
            print('File selected:   ', input_spec_list, '(only 1 file).',
                  'You set z='+"{:.4f}".format(z_guess),'for this object')
        
    try:
        dic_templates[templatename]
        print(    'Template:        ', templatename, 'is selected.')
    except KeyError:
        templatename = 'vvds_spiral'
        print(    '[Error] Template:', templatename, 'does not exist. Switched to default template: "VVDS Spiral".')   
else: # No input arguments
    spec1dfolderpath = 'input_spec1d/'
    spec1dobjname    = 'spec1d.m46.071.A2552.fits'
    spec1dfilepath   = spec1dfolderpath + spec1dobjname
    input_spec_list  = [spec1dfilepath]
    templatename     = 'vvds_spiral'
    z_guess          = 0.5135
    print('No input arguments detected. Begin with default:', spec1dfilepath)
    
os.makedirs('redshifts_JPGs', exist_ok=True)


# FOR EVERY SPEC1D OBJECT...
def main_func(input_spec_list, spec1dfolderpath, templatename, z_guess):
    for spec1dfilepath in input_spec_list:
        spec1dobjname = spec1dfilepath[len(spec1dfolderpath):]
        
        # **observed** spectrum (obs_spec)
        hdulist = fits.open(spec1dfilepath)
        print('Opened...         ' + spec1dfilepath)
        hdr01data = hdulist[1].data
        w = hdr01data["LAMBDA"][0]
        f = hdr01data["FLUX"][0]
        e = hdr01data["IVAR"][0]
        arr_obs_wave = binning_1d(w, binfactor) * u.AA
        arr_obs_flux = binning_1d(f, binfactor) * u.dimensionless_unscaled
        arr_obs_erro = binning_1d(e, binfactor) * u.dimensionless_unscaled
        del e, f, w
        
        # **template** spectrum (tem_spec)
        arr_tem_wave = dic_templates[templatename][0] * u.AA
        arr_tem_flux = dic_templates[templatename][1] * u.dimensionless_unscaled
        arr_tem_erro = 0.2e-19 * np.ones(len(arr_tem_flux)) * u.dimensionless_unscaled
        
        
        """
        FLUX-CONTINUUM NORMALIZATION
        
        # 1) Estimate continuum by ignoring all obscured lines
             "//2*2 + 1": kernal_size must be odd
                     "5": self-def. Smaller -> fewer obscured lines
                   
        # 2) Calculate relative (template-continuum - obs-flux-continuum)
        
        # 3) Fit continuum -
             Assemble [wavelength, spec, error] into the spec objects
             spectral_axis: array of float, quantity with unit u.AA
             flux:          array of float, quantity with unit u.Jy, 
                            u.Unit("erg / (Angstrom s cm2)"), or defined
             uncertainty:   array of float, quantity with unit u.Jy, 
                            u.Unit("erg / (Angstrom s cm2)"), or defined
             This step gives you a redshift probability distribution. 
             Only redshifts > 0 are realistic.
             
        # 4) Amplitude normalization
             
        """
        arr_obs_cont_flux = signal.medfilt(arr_obs_flux, 
                                           kernel_size=int((len(arr_obs_flux)/5)//2*2 + 1)
                                           )
        arr_tem_cont_flux = signal.medfilt(arr_tem_flux, 
                                           kernel_size=int((len(arr_obs_flux)/5)//2*2 + 1)
                                           )
        arr_obs_rela_flux = arr_obs_flux - arr_obs_cont_flux
        arr_tem_rela_flux = arr_tem_flux - arr_tem_cont_flux
        obs_rela_spec = Spectrum1D(spectral_axis = arr_obs_wave, 
                                   flux          = arr_obs_rela_flux, 
                                   uncertainty   = StdDevUncertainty(arr_obs_erro))
        tem_rela_spec = Spectrum1D(spectral_axis = arr_tem_wave, 
                                   flux          = arr_tem_rela_flux, 
                                   uncertainty   = StdDevUncertainty(arr_tem_erro))
        corr, lag = correlation.template_correlate(observed_spectrum = obs_rela_spec, 
                                                   template_spectrum = tem_rela_spec,
                                                   lag_units=u.dimensionless_unscaled)
        corr = corr[lag > 0.]
        lag  =  lag[lag > 0.]
        corr_normalized = corr/sum(corr)
        z_peak = lag[np.argmax(corr_normalized)].value
        if z_guess == None:
            z_guess = z_peak
        
        normal_fac = sum(arr_obs_rela_flux) / sum(arr_tem_rela_flux)
        arr_tem_rela_flux *= normal_fac
        arr_tem_wave      *= (1+z_guess)
        
        # Resamples 1): We want a new obs_flux taken from tem_wave's axis values
        # Resamples 2): make f_obs(x_obs)
        interpo_func_obs_cont = interpolate.interp1d(arr_obs_wave, arr_obs_cont_flux, kind='linear')
        # Resamples 3): make common wavelength axis w/ same numbers
        mask = (arr_tem_wave >= min(arr_obs_wave)) & (arr_tem_wave <= max(arr_obs_wave))
        arr_temobs_common_wave   =      arr_tem_wave[mask]
        arr_tem_cont_flux_new    = arr_tem_cont_flux[mask]
        # Resamples 4): apply f_obs(x) for x in x_tem
        arr_obs_cont_flux_new    = interpo_func_obs_cont(arr_temobs_common_wave)
        
        arr_obs_tem_cont_ratio  = arr_obs_cont_flux_new / arr_tem_cont_flux_new
        arr_tem_cont_flux_new   = arr_tem_cont_flux[mask] * arr_obs_tem_cont_ratio
        arr_tem_flux_new        = arr_tem_flux[mask] * arr_obs_tem_cont_ratio
        arr_tem_wave_new        = arr_tem_wave[mask]
        
        
        
        # All plottings
        fig = plt.figure(figsize=(12,12),dpi=200)
        plt.subplots_adjust(hspace=0.35, wspace=0.2) # h=height
        gs = fig.add_gridspec(3, 2,
                              height_ratios=[3,2,2],
                              width_ratios=[1, 1])
        ax1 = fig.add_subplot(gs[0, :])
        ax4 = fig.add_subplot(gs[1, :])
        ax2 = fig.add_subplot(gs[2, 0])
        ax3 = fig.add_subplot(gs[2, 1])
        if z_guess != z_peak:
            ax1.set_title(spec1dobjname+
                      ' --- z='+"{:.4f}".format(z_guess)+
                      ' (manual)'
                      , fontsize=20, loc='left')
        else:
            ax1.set_title(spec1dobjname+
                      ' --- z='+"{:.4f}".format(z_guess)
                      , fontsize=20, loc='left')
        
        # Sub-plot for observed and template spectra
        # ax1.grid(alpha=0.75, linestyle='--', zorder=-1)
        ax1.plot(arr_obs_wave.value,
                 arr_obs_flux.value, 
                 label='Observed spectrum',
                 color='black', linewidth=1, zorder=11)
        ax1.plot(arr_tem_wave_new.value, 
                 arr_tem_flux_new.value, 
                 label='Template: '+templatename+' at z='+"{:.4f}".format(z_guess),
                 color='red', linewidth=1,  zorder=11, linestyle='-')
        ax1.plot(arr_obs_wave.value, 
                 arr_obs_cont_flux, 
                 label='Observed continuum', alpha=0.75,
                 color='green', linewidth=2, zorder=10, linestyle=':')
        ax1.plot(arr_tem_wave_new.value, 
                 arr_tem_cont_flux_new, 
                 label='Template continuum', alpha=0.75, 
                 color='red', linewidth=2, zorder=10, linestyle=':')
        ax1.fill_between(arr_obs_wave.value,
                         arr_obs_erro.value,
                         np.zeros(len(arr_obs_flux)), 
                         color='black', linewidth=0, alpha=.15, zorder=9)
        if templatename in ['vvds_elliptical','vvds_s0','red_galaxy']:
            alpha_abs = 0.75
            alpha_emi = 0.5
        else:
            alpha_abs = 0.5
            alpha_emi = 0.75
        xmin1 = arr_obs_wave[ 0].value
        xmax1 = arr_obs_wave[-1].value
        ymin1 = ax1.get_ylim()[0]
        ymax1 = ax1.get_ylim()[1]
        plot_lines(dic_emi_lines, dic_abs_lines, z_guess, ax1, 
                   alpha_abs, alpha_emi, xmin1, xmax1, ymin1, ymax1)
        ax1.set_ylabel('Flux')
        ax1.minorticks_on()
        ax1.tick_params(  which='minor', length=3,bottom=True,top=True,direction='in')
        ax1.tick_params(  which='major', length=6,bottom=True,top=True,direction='in')
        ax1.set_xlabel(r'Observed Wavelength ($\AA$)')
        ax1.set_xlim(xmin1, xmax1)
        ax1.legend(loc='upper left', prop={'size': 10})
        
        # Sub-plot for subtracted spectra
        # ax4.grid(alpha=0.75, linestyle='--', zorder=-1)
        ax4.plot(arr_obs_wave,
                 arr_obs_rela_flux, 
                 label='Observed', 
                 color='black', linewidth=1, zorder=10)
        ax4.plot(arr_tem_wave, 
                 arr_tem_rela_flux, 
                 label='Template', 
                 color='red', linewidth=1, zorder=9, linestyle='-')
        ymin4 = ax4.get_ylim()[0]
        ymax4 = ax4.get_ylim()[1]
        plot_lines(dic_emi_lines, dic_abs_lines, z_guess, ax4, 
                   alpha_abs, alpha_emi, xmin1, xmax1, ymin4, ymax4)
        ax4.set_ylabel('Relative Flux (spectrum - continuum)')
        ax4.minorticks_on()
        ax4.tick_params(which='minor', length=3,bottom=True,top=True,direction='in')
        ax4.tick_params(which='major', length=6,bottom=True,top=True,direction='in',
                        axis='x', labeltop=True, labelbottom=True)
        ax4.set_xlabel(r'Observed Wavelength ($\AA$)')
        ax4.set_xlim(xmin1, xmax1)
        ax4.legend(loc='upper left', prop={'size': 10})
        
        # Sub-plot for redshift probability distribution and guess
        ax2.grid(alpha=0.75, linestyle='--', zorder=-1)
        ax2.plot(lag, 
                corr_normalized, 
                color='green', linewidth=1, zorder=10)
        ax2.set_xlim(left=0, right=None)
        ax2.set_ylim(bottom=0, top=None)
        
        # Sub-plot for zoom-in redshift peak
        try:
            ax3.grid(alpha=0.5, linestyle='--', zorder=-1)
            ax3.plot(lag, 
                    corr_normalized, 
                    color='green', linewidth=1, zorder=10)
            z_peak = lag[np.argmax(corr_normalized)]
            ymin  = np.min(corr_normalized)
            ymax  = np.max(corr_normalized)*1.25
            xmin2 = np.min(lag[corr_normalized > np.max(corr_normalized)*0.5]) - 0.02
            xmax2 = np.max(lag[corr_normalized > np.max(corr_normalized)*0.5]) + 0.02
            ymin2 = np.max(corr_normalized)*0.0
            ymax2 = np.max(corr_normalized)*1.25
            ax3.vlines(x=z_peak,  ymin=ymin2, ymax=ymax2, 
                       ls='--', zorder=8, color='black', linewidth=1)
            ax3.text(z_peak, ymax2, r"$z_\mathrm{peak}=$"+"{:.4f}".format(z_peak), size=10, 
                     horizontalalignment='center', verticalalignment='bottom', 
                     color='green', zorder=7)
            if z_guess != z_peak:
                ax2.vlines(x=z_guess, ymin=ymin, ymax=ymax, 
                           ls='-',  zorder=6, color='red', linewidth=1)
                ax3.vlines(x=z_guess, ymin=ymin2, ymax=ymin2+(ymax2-ymin2)*0.9, 
                           ls='-',  zorder=6, color='red', linewidth=1)
                ax3.text(z_guess, ymin2+(ymax2-ymin2)*0.98,  "{:.4f}".format(z_guess), size=10, 
                         horizontalalignment='center', verticalalignment='top', 
                         color='red', zorder=7)
            else:
                ax3.vlines(x=z_peak, ymin=ymin2, ymax=ymax2, 
                           ls='-',  zorder=8, color='red', alpha=0.5, linewidth=1)
            ax3.set_xlim(xmin2, xmax2)
            ax3.set_ylim(ymin2, ymax2)
        except TypeError:
            ax3.text(0.5, 0.5, 'Zoom-in Unavailable', size=20,
                     horizontalalignment='center', verticalalignment='center', 
                     color='red', zorder=7,
                     transform=ax3.transAxes)
            print('[Error] Zoom-in failed: check if any index is invalid!')
            pass
        ax2.set_ylabel('Correlation')
        ax3.set_ylabel('Correlation')
        ax2.minorticks_on()
        ax3.minorticks_on()
        ax2.tick_params(which='minor', length=3, bottom=True,top=True, direction='in')
        ax2.tick_params(which='major', length=9, bottom=True,top=True, direction='in',
                        axis='x', labeltop='on', labelbottom='on')
        ax3.tick_params(which='minor', length=3, bottom=True,top=True, direction='in')
        ax3.tick_params(which='major', length=15,bottom=True,top=True, direction='inout',
                        axis='x', labeltop='on', labelbottom='on', labelsize=12)
        ax2.set_xlabel(r'Redshift ($z$)')
        ax3.set_xlabel(r'Redshift ($z$)')
        # ax3.set_yscale('log')
        if z_guess != z_peak:
            plt.savefig("redshifts_JPGs/"+spec1dobjname+"_adjusted_new.jpg", dpi=300, bbox_inches='tight')
        else:
            plt.savefig("redshifts_JPGs/"+spec1dobjname+".jpg", dpi=300, bbox_inches='tight')
        z_guess = None
    return print('Complete:         '+spec1dfilepath)

main_func(input_spec_list, spec1dfolderpath, templatename, z_guess)