import numpy as np
import sys
import os
import matplotlib.pyplot as plt

from codes.read_spec1d        import read_spec1d
from codes.read_spec2d        import read_spec2d
from codes.plots_spec_z_combo import plots_spec_z_combo
from codes.fitting            import cross_corr
from codes.binning_1d         import binning_1d

from codes.load_templates     import read_templates
from codes.load_lines         import read_lines
global dic_templates, dic_emi_lines, dic_abs_lines
dic_templates                = read_templates()
dic_emi_lines, dic_abs_lines = read_lines()
print('Template and lines loading complete.')

z_guess1 = None
z_guess2 = None

binfactor_spec = 30
binfactor_temp = 5

# Turn off: ``RuntimeWarning: invalid value encountered in divide.''
np.seterr(divide='ignore', invalid='ignore')

# USER INPUTS...
# $ python firefly.py input_spec1d/ input_spec2d/ output_JPGs/ -all sdss_luminous_red sdss_late_type
# $ python firefly.py input_spec1d/ input_spec2d/ output_JPGs/ spec1d.m.030.A.fits sdss_luminous_red sdss_late_type 0.68 0.67
inputlist = sys.argv

if len(inputlist) > 1:
    spec1dfolderpath = inputlist[1]
    spec2dfolderpath = inputlist[2]
    outputfolderpath = inputlist[3]
    spec1dobjnames   = inputlist[4]
    templatename1    = inputlist[5]
    templatename2    = inputlist[6]
    # check illegal inputs
    if spec1dfolderpath[-1] != '/':
        spec1dfolderpath += '/'
    print('Spec1d Folder:   ', spec1dfolderpath, 'is specified.')
    if spec2dfolderpath[-1] != '/':
        spec2dfolderpath += '/'
    print('Spec2d Folder:   ', spec2dfolderpath, 'is specified.')
    if outputfolderpath[-1] != '/':
        outputfolderpath += '/'
    print('Output Folder:   ', outputfolderpath, 'is specified.')
    #
    if spec1dobjnames == '-all':
        input_spec1d_lst = []
        short_list      = os.listdir(spec1dfolderpath)
        for f in short_list: input_spec1d_lst += [spec1dfolderpath + f]
        print('File selected:   ', len(input_spec1d_lst), 'files in this spec1d folder.')
    else:
        input_spec1d_lst = [spec1dfolderpath + spec1dobjnames]
        z_guess1 = float(inputlist[7])
        z_guess2 = float(inputlist[8])
        print('File selected:   ', input_spec1d_lst, '(only 1 file).',
              'You set z='+"{:.4f} {:.4f}".format(z_guess1,z_guess2),'for this object')
    #
    try:
        dic_templates[templatename1]
        print('Template 1:      ', templatename1, 'is selected.')
    except KeyError:
        templatename1 = 'vvds_spiral'
        print('[Error] Template:', templatename1, 'does not exist.',
              'Available templates: ', list(dic_templates.keys()),
              'Switched to default template: "VVDS Spiral".')   
    #
    try:
        dic_templates[templatename2]
        print('Template 2:      ', templatename2, 'is selected.')
    except KeyError:
        templatename2 = 'vvds_elliptical'
        print('[Error] Template:', templatename2, 'does not exist.',
              'Available templates: ', list(dic_templates.keys()),
              'Switched to default template: "VVDS Elliptical".')
else: # No input arguments
    spec1dfolderpath = 'input_spec1d_TEST/'
    input_spec1d_lst = []
    short_list       = os.listdir(spec1dfolderpath)
    for f in short_list: input_spec1d_lst += [spec1dfolderpath + f]
    spec2dfolderpath = 'input_spec2d_TEST/'
    # spec2dobjname    = 'spec2d.m46.077.A2552.fits'
    # spec2dfilepath   = spec2dfolderpath + spec2dobjname
    # input_spec2d_lst = [spec2dfilepath]
    templatename1    = 'sdss_luminous_red'
    templatename2    = 'sdss_late_type'
    z_guess1         = None
    z_guess2         = None
    
os.makedirs('redshifts_JPGs', exist_ok=True)


# FOR EVERY SPEC1D OBJECT...
def main_func(input_spec1d_lst, spec1dfolderpath, spec2dfolderpath,
              templatename1, templatename2, z_guess1, z_guess2,
              binfactor_spec, binfactor_temp):
    
    for spec1dfilepath in input_spec1d_lst:
        spec1dobjname = spec1dfilepath[len(spec1dfolderpath):]
        input_spec2d_lst = []
        short_list       = os.listdir(spec2dfolderpath)
        for f in short_list: input_spec2d_lst += [spec2dfolderpath + f]
        try: # slit number in spec2d list
            slitnumber = spec1dobjname.split('.')[2]
            for spec2dfilepath in input_spec2d_lst:
                spec2dobjname = spec2dfilepath[len(spec2dfolderpath):]
                if spec2dobjname.split('.')[2] == slitnumber:
                    break
        except:
            print('Spec-2d not found. Check if any filename does not have typos.')
        
        # 2d spectrum
        farr,vmin,vmax,warr = read_spec2d(spec2dfilepath)
        
        # **observed** spectrum --------------------------------------------------------------
        arr_obs_wave, arr_obs_flux, arr_obs_erro = read_spec1d(spec1dfilepath, binfactor_spec)
        
        wave_min, wave_max = np.min(arr_obs_wave), np.max(arr_obs_wave)
        wave_range         = wave_max - wave_min
        
        mask_positive_flux = arr_obs_flux > 0
        mask_left_edge     = arr_obs_wave > (wave_min + 0.03*wave_range)
        mask_right_edge    = arr_obs_wave < (wave_min + 0.97*wave_range)
        mask_atm_absorp_Al = arr_obs_wave < 7586
        mask_atm_absorp_Ar = arr_obs_wave > 7708
        mask_atm_absorp_Bl = arr_obs_wave < 6864
        mask_atm_absorp_Br = arr_obs_wave > 6945
        mask_all_1 = mask_positive_flux & mask_left_edge & mask_right_edge
        mask_all_2 = mask_atm_absorp_Al | mask_atm_absorp_Ar
        mask_all_3 = mask_atm_absorp_Bl | mask_atm_absorp_Br
        mask_all   = mask_all_1 & mask_all_2 & mask_all_3
        
        arr_obs_flux = arr_obs_flux[mask_all]
        arr_obs_wave = arr_obs_wave[mask_all]
        arr_obs_erro = arr_obs_erro[mask_all]
        
        arr_obs = np.concatenate(([arr_obs_wave],
                                  [arr_obs_flux],
                                  [arr_obs_erro]), axis=0)
        
        # **template** spectrum (tem_spec) ----------------------------------------------------------------------------
        # No.1
        arr_tem1wave = binning_1d(dic_templates[templatename1][0], binfactor_temp)
        arr_tem1flux = binning_1d(dic_templates[templatename1][1], binfactor_temp)
        if len(dic_templates[templatename1]) <= 2:
            arr_tem1erro = binning_1d(0.2e-19 * np.ones(len(arr_tem1flux)), binfactor_temp)
        else:
            arr_tem1erro = binning_1d(dic_templates[templatename1][2], binfactor_temp)
        arr_tem1 = np.concatenate(([arr_tem1wave],
                                   [arr_tem1flux],
                                   [arr_tem1erro]), axis=0)
        # No.2
        arr_tem2wave = binning_1d(dic_templates[templatename2][0], binfactor_temp)
        arr_tem2flux = binning_1d(dic_templates[templatename2][1], binfactor_temp)
        if len(dic_templates[templatename2]) <= 2:
            arr_tem2erro = binning_1d(0.2e-19 * np.ones(len(arr_tem2flux)), binfactor_temp)
        else:
            arr_tem2erro = binning_1d(dic_templates[templatename2][2], binfactor_temp)
        arr_tem2 = np.concatenate(([arr_tem2wave],
                                   [arr_tem2flux],
                                   [arr_tem2erro]), axis=0)
        
        # Fittings -----------------------------------------------------------------------
        (z_peak1, z_guess1,
        lag1, corr1_normalized, 
        arr_obs_cont,        # Figure 1+3: 2/4 (1/4=arr_obs)
        arr_tem1_masked,     # Figure 1+3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
        arr_tem1cont_masked, # Figure 1+3, 4/4
        arr_obs_rela,    # Figure 2, 1/2
        arr_tem1rela,    # Figure 2, 2/2
        arr_tem1_masked_newz,     # Figure 3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
        arr_tem1cont_masked_newz, # Figure 3, 4/4
        arr_tem1relaC             # Figure 4, 2/2
        ) = cross_corr(arr_obs, arr_tem1, 
                       z_guess1, dic_emi_lines)
        
        (z_peak2, z_guess2,
        lag2, corr2_normalized, 
        arr_obs_cont,        # Figure 1+3: 2/4 (1/4=arr_obs)
        arr_tem2_masked,     # Figure 1+3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
        arr_tem2cont_masked, # Figure 1+3, 4/4
        arr_obs_rela,    # Figure 2, 1/2
        arr_tem2rela,    # Figure 2, 2/2
        arr_tem2_masked_newz,     # Figure 3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
        arr_tem2cont_masked_newz, # Figure 3, 4/4
        arr_tem2relaC             # Figure 4, 2/2
        ) = cross_corr(arr_obs, arr_tem2, 
                       z_guess2, dic_emi_lines)
        
        # All plottings
        fig = plt.figure(figsize=(22,12),dpi=100)
        spec2dbox_ratio = (fig.get_size_inches()[0]/2 / (fig.get_size_inches()[1]/(11)))
        plt.subplots_adjust(hspace=0.15, wspace=0.20) # h=height
        gs = fig.add_gridspec(6, 4,
                              height_ratios=[1,2,2,2,2,2],
                              width_ratios=[1,1,1,1])
        ax_2dspecimg_lf = fig.add_subplot(gs[0, 0:2])
        ax_orig_spec_lf = fig.add_subplot(gs[1, 0:2])
        ax_rela_spec_lf = fig.add_subplot(gs[2, 0:2])
        ax_orig_spe2_lf = fig.add_subplot(gs[3, 0:2])
        ax_rela_spe2_lf = fig.add_subplot(gs[4, 0:2])
        ax_corr_dist_lf = fig.add_subplot(gs[5, 0])
        ax_peak_dist_lf = fig.add_subplot(gs[5, 1])
        
        plots_spec_z_combo( ax_2dspecimg_lf,
                            farr, vmin, vmax, warr, spec2dbox_ratio,
                            ax_orig_spec_lf, ax_rela_spec_lf, 
                            ax_corr_dist_lf, ax_peak_dist_lf,
                            ax_orig_spe2_lf, ax_rela_spe2_lf,
                            dic_abs_lines, dic_emi_lines,
                            spec1dobjname, templatename1, z_peak1, z_guess1,
                            lag1, corr1_normalized, 
                            arr_obs,             # Figure 1+3: 1/4=arr_obs
                            arr_obs_cont,        # Figure 1+3: 2/4 (1/4=arr_obs)
                            arr_tem1_masked,     # Figure 1: 3/4 (Can't use arr_tem1 b/c un-normalized!)
                            arr_tem1cont_masked, # Figure 1, 4/4
                            arr_obs_rela,    # Figure 2+4, 1/2
                            arr_tem1rela,    # Figure 2, 2/2
                            arr_tem1_masked_newz,     # Figure 3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
                            arr_tem1cont_masked_newz, # Figure 3, 4/4
                            arr_tem1relaC,            # Figure 4, 2/2
                            ['red', 'red'])
        
        ax_2dspecimg_rt = fig.add_subplot(gs[0, 2:])
        ax_orig_spec_rt = fig.add_subplot(gs[1, 2:])
        ax_rela_spec_rt = fig.add_subplot(gs[2, 2:])
        ax_orig_spe2_rt = fig.add_subplot(gs[3, 2:])
        ax_rela_spe2_rt = fig.add_subplot(gs[4, 2:])
        ax_corr_dist_rt = fig.add_subplot(gs[5, 2])
        ax_peak_dist_rt = fig.add_subplot(gs[5, 3])
        
        plots_spec_z_combo( ax_2dspecimg_rt,
                            farr, vmin, vmax, warr, spec2dbox_ratio,
                            ax_orig_spec_rt, ax_rela_spec_rt, 
                            ax_corr_dist_rt, ax_peak_dist_rt,
                            ax_orig_spe2_rt, ax_rela_spe2_rt,
                            dic_abs_lines, dic_emi_lines,
                            spec1dobjname, templatename2, z_peak2, z_guess2,
                            lag2, corr2_normalized, 
                            arr_obs,             # Figure 1+3: 1/4=arr_obs
                            arr_obs_cont,        # Figure 1+3: 2/4 (1/4=arr_obs)
                            arr_tem2_masked,     # Figure 1: 3/4 (Can't use arr_tem1 b/c un-normalized!)
                            arr_tem2cont_masked, # Figure 1, 4/4
                            arr_obs_rela,    # Figure 2+4, 1/2
                            arr_tem2rela,    # Figure 2, 2/2
                            arr_tem2_masked_newz,     # Figure 3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
                            arr_tem2cont_masked_newz, # Figure 3, 4/4
                            arr_tem2relaC,            # Figure 4, 2/2
                            ['magenta', 'magenta'])
        
        
        fig.suptitle("Slit {} ({} + {})".format(slitnumber,spec1dobjname, spec2dobjname), 
                     x=0.5, y=0.92, fontsize=18, va='top')
        ax_2dspecimg_lf.set_title('Template 1: '+templatename1,
                                  fontsize=14, loc='center')
        ax_2dspecimg_rt.set_title('Template 2: '+templatename2,
                                  fontsize=14, loc='center')
        
        if z_guess1 != z_peak1 or z_guess2 != z_peak2:
            plt.savefig(outputfolderpath+spec1dobjname+"_adjusted_test.jpg", dpi=150, bbox_inches='tight')
        else:
            plt.savefig(outputfolderpath+spec1dobjname+".jpg", dpi=150, bbox_inches='tight')
        z_guess1 = None
        z_guess2 = None
        
        print('Complete:         '+spec1dfilepath)
        
    return 

main_func(input_spec1d_lst, spec1dfolderpath, spec2dfolderpath, 
          templatename1, templatename2, z_guess1, z_guess2,
          binfactor_spec, binfactor_temp)