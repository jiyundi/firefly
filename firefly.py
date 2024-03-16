import numpy as np
from astropy import units as u
import sys
import os
import matplotlib.pyplot as plt

from code.read_spec1d        import read_spec1d
from code.plots_spec_z_combo import plots_spec_z_combo
from code.fitting            import cross_corr

from code.load_templates     import read_templates
from code.load_lines         import read_lines
global dic_templates, dic_emi_lines, dic_abs_lines
dic_templates                = read_templates()
dic_emi_lines, dic_abs_lines = read_lines()
print('Template and lines loading complete.')

z_guess1 = None
z_guess2 = None

# Turn off: ``RuntimeWarning: invalid value encountered in divide.''
np.seterr(divide='ignore', invalid='ignore')

# USER INPUTS...
inputlist = sys.argv

if len(inputlist) > 1:
    spec1dfolderpath = inputlist[1]
    spec1dobjnames   = inputlist[2]
    templatename1    = inputlist[3]
    templatename2    = inputlist[4]
    # check illegal inputs
    if spec1dfolderpath[-1] != '/':
        spec1dfolderpath += '/'
        print('Spec1d Folder:   ', spec1dfolderpath, 'is specified.')
    try:
        dic_templates[templatename1]
        print('Template 1:      ', templatename1, 'is selected.')
    except KeyError:
        templatename1 = 'vvds_spiral'
        print('[Error] Template:', templatename1, 'does not exist.',
              'Available templates: ', list(dic_templates.keys()),
              'Switched to default template: "VVDS Spiral".')   
    try:
        dic_templates[templatename2]
        print('Template 2:      ', templatename2, 'is selected.')
    except KeyError:
        templatename2 = 'vvds_elliptical'
        print('[Error] Template:', templatename2, 'does not exist.',
              'Available templates: ', list(dic_templates.keys()),
              'Switched to default template: "VVDS Elliptical".')
    if spec1dobjnames == '-all':
        input_spec_list = []
        short_list      = os.listdir(spec1dfolderpath)
        for f in short_list: input_spec_list += [spec1dfolderpath + f]
        print('File selected:   ', len(input_spec_list), 'files in this spec1d folder.')
    else:
        input_spec_list = [spec1dfolderpath + spec1dobjnames]
        z_guess1 = float(inputlist[5])
        z_guess2 = float(inputlist[6])
        print('File selected:   ', input_spec_list, '(only 1 file).',
              'You set z='+"{:.4f} {:.4f}".format(z_guess1,z_guess2),'for this object')
else: # No input arguments
    spec1dfolderpath = 'input_spec1d/'
    spec1dobjname    = 'spec1d.m46.077.A2552.fits'
    spec1dfilepath   = spec1dfolderpath + spec1dobjname
    input_spec_list  = [spec1dfilepath]
    templatename1    = 'vvds_spiral'
    templatename2    = 'vvds_elliptical'
    z_guess1         = 0.3968
    z_guess2         = 0.8590
    print('No input arguments detected. Begin with default:', spec1dfilepath)
    
os.makedirs('redshifts_JPGs', exist_ok=True)


# FOR EVERY SPEC1D OBJECT...
def main_func(input_spec_list, spec1dfolderpath, 
              templatename1, templatename2, z_guess1, z_guess2):
    
    for spec1dfilepath in input_spec_list:
        spec1dobjname = spec1dfilepath[len(spec1dfolderpath):]
        
        # **observed** spectrum
        arr_obs_wave, arr_obs_flux, arr_obs_erro = read_spec1d(spec1dfilepath)
        
        arr_obs_flux_raw = arr_obs_flux.copy()
        arr_obs_wave_raw = arr_obs_wave.copy()
        arr_obs_erro_raw = arr_obs_erro.copy()
        
        wave_min, wave_max = np.min(arr_obs_wave), np.max(arr_obs_wave)
        wave_range         = wave_max - wave_min
        
        mask_positive_flux = arr_obs_flux > 0
        mask_left_edge     = arr_obs_wave > (wave_min + 0.03*wave_range)
        mask_right_edge    = arr_obs_wave < (wave_min + 0.97*wave_range)
        mask_atm_absorp_Al = arr_obs_wave.value < 7586
        mask_atm_absorp_Ar = arr_obs_wave.value > 7708
        mask_atm_absorp_Bl = arr_obs_wave.value < 6864
        mask_atm_absorp_Br = arr_obs_wave.value > 6945
        mask_all_1 = mask_positive_flux & mask_left_edge & mask_right_edge
        mask_all_2 = mask_atm_absorp_Al | mask_atm_absorp_Ar
        mask_all_3 = mask_atm_absorp_Bl | mask_atm_absorp_Br
        mask_all   = mask_all_1 & mask_all_2 & mask_all_3
        
        arr_obs_flux = arr_obs_flux[mask_all]
        arr_obs_wave = arr_obs_wave[mask_all]
        arr_obs_erro = arr_obs_erro[mask_all]
        
        # **template** spectrum (tem_spec)
        arr_tem1wave = dic_templates[templatename1][0] * u.AA
        arr_tem1flux = dic_templates[templatename1][1] * u.dimensionless_unscaled
        arr_tem1erro = 0.2e-19 * np.ones(len(arr_tem1flux)) * u.dimensionless_unscaled
        
        arr_tem2wave = dic_templates[templatename2][0] * u.AA
        arr_tem2flux = dic_templates[templatename2][1] * u.dimensionless_unscaled
        arr_tem2erro = 0.2e-19 * np.ones(len(arr_tem2flux)) * u.dimensionless_unscaled
        
        # FLUX-CONTINUUM NORMALIZATION AND SUBSTRACTION
        result_dic = cross_corr(arr_obs_wave, arr_obs_flux, arr_obs_erro,
                                arr_tem1wave, arr_tem1flux, arr_tem1erro,
                                z_guess1)
        z_guess1              = result_dic['z_guess']
        z_peak1               = result_dic['z_peak']
        lag1                  = result_dic['lag']
        corr_normalized1      = result_dic['corr_normalized']
        arr_obs_rela_flux     = result_dic['arr_obs_rela_flux']
        arr_obs_cont_wave     = result_dic['arr_obs_cont_wave']
        arr_obs_cont_flux     = result_dic['arr_obs_cont_flux']
        arr_tem1rela_flux     = result_dic['arr_tem_rela_flux']
        arr_tem1wave_new      = result_dic['arr_tem_wave_new']
        arr_tem1flux_new      = result_dic['arr_tem_flux_new']
        arr_tem1cont_flux_new = result_dic['arr_tem_cont_flux_new']
        
        result_dic = cross_corr(arr_obs_wave, arr_obs_flux, arr_obs_erro,
                                arr_tem2wave, arr_tem2flux, arr_tem2erro,
                                z_guess2)
        z_guess2              = result_dic['z_guess']
        z_peak2               = result_dic['z_peak']
        lag2                  = result_dic['lag']
        corr_normalized2      = result_dic['corr_normalized']
        arr_obs_rela_flux     = result_dic['arr_obs_rela_flux']
        arr_obs_cont_wave     = result_dic['arr_obs_cont_wave']
        arr_obs_cont_flux     = result_dic['arr_obs_cont_flux']
        arr_tem2rela_flux     = result_dic['arr_tem_rela_flux']
        arr_tem2wave_new      = result_dic['arr_tem_wave_new']
        arr_tem2flux_new      = result_dic['arr_tem_flux_new']
        arr_tem2cont_flux_new = result_dic['arr_tem_cont_flux_new']
        
        # All plottings
        fig = plt.figure(figsize=(16,10),dpi=150)
        plt.subplots_adjust(hspace=0.4, wspace=0.35) # h=height
        gs = fig.add_gridspec(3, 4,
                              height_ratios=[3,2,2],
                              width_ratios=[1,1,1,1])
        ax1 = fig.add_subplot(gs[0, 0:2])
        ax4 = fig.add_subplot(gs[1, 0:2])
        ax2 = fig.add_subplot(gs[2, 0])
        ax3 = fig.add_subplot(gs[2, 1])
        
        plots_spec_z_combo( ax1, ax2, ax3, ax4,
                            dic_abs_lines, dic_emi_lines,
                            spec1dobjname, templatename1, z_peak1, z_guess1,
                            lag1, corr_normalized1,
                            arr_obs_wave_raw,  arr_obs_flux_raw,  arr_obs_erro_raw,
                            arr_obs_wave,      arr_obs_rela_flux,
                            arr_obs_cont_wave, arr_obs_cont_flux, 
                            arr_tem1wave,      arr_tem1rela_flux,
                            arr_tem1wave_new,  arr_tem1flux_new,  arr_tem1cont_flux_new)
        
        ax5 = fig.add_subplot(gs[0, 2:])
        ax8 = fig.add_subplot(gs[1, 2:])
        ax6 = fig.add_subplot(gs[2, 2])
        ax7 = fig.add_subplot(gs[2, 3])
        
        plots_spec_z_combo( ax5, ax6, ax7, ax8,
                            dic_abs_lines, dic_emi_lines,
                            spec1dobjname, templatename2, z_peak2, z_guess2,
                            lag2, corr_normalized2,
                            arr_obs_wave_raw,  arr_obs_flux_raw,  arr_obs_erro_raw,
                            arr_obs_wave,      arr_obs_rela_flux,
                            arr_obs_cont_wave, arr_obs_cont_flux, 
                            arr_tem2wave,      arr_tem2rela_flux,
                            arr_tem2wave_new,  arr_tem2flux_new,  arr_tem2cont_flux_new)
        
        fig.suptitle(spec1dobjname, x=0.5, y=0.93, fontsize=18, va='top')
        ax1.set_title('Template 1: '+templatename1, 
                      fontsize=14, loc='center')
        ax5.set_title('Template 2: '+templatename2, 
                      fontsize=14, loc='center')
        
        if z_guess1 != z_peak1 or z_guess2 != z_peak2:
            plt.savefig("redshifts_JPGs/"+spec1dobjname+"_adjusted_test.jpg", dpi=150, bbox_inches='tight')
        else:
            plt.savefig("redshifts_JPGs/"+spec1dobjname+".jpg", dpi=150, bbox_inches='tight')
        z_guess1 = None
        z_guess2 = None
        
        print('Complete:         '+spec1dfilepath)
        
    return 

main_func(input_spec_list, spec1dfolderpath, 
          templatename1, templatename2, z_guess1, z_guess2)