import numpy as np
from codes.plot_lines     import plot_lines
from codes.plot_blocks    import plot_blocks, plot_nega_flux, plot_nega_flux_rela

def plots_spec_z_combo(ax1, ax2, ax3, ax4,
                       dic_abs_lines, dic_emi_lines,
                       spec1dobjname, templatename, z_peak, z_guess,
                       lag, corr_normalized,
                       arr_obs_wave_raw,  arr_obs_flux_raw,  arr_obs_erro_raw,
                       arr_obs_wave,      arr_obs_rela_flux,
                       arr_obs_cont_wave, arr_obs_cont_flux, 
                       arr_tem_wave,      arr_tem_rela_flux,
                       arr_tem_wave_new,  arr_tem_flux_new,  arr_tem_cont_flux_new):
    # Sub-plot for observed and template spectra
    ax1.plot(arr_obs_wave_raw.value,
             arr_obs_flux_raw.value, 
             label='Observed spectrum',
             color='black', linewidth=1, zorder=11)
    ax1.plot(arr_obs_cont_wave.value, 
             arr_obs_cont_flux, 
             label='Observed continuum', alpha=0.75,
             color='green', linewidth=2, zorder=10, linestyle=':')
    ax1.plot(arr_tem_wave_new.value, 
             arr_tem_flux_new.value, 
             label='Template at z='+"{:.4f}".format(z_guess),
             color='red', linewidth=1,  zorder=11, linestyle='-')
    ax1.plot(arr_tem_wave_new.value, 
             arr_tem_cont_flux_new, 
             label='Template continuum', alpha=0.75, 
             color='red', linewidth=2, zorder=10, linestyle=':')
    ax1.fill_between(arr_obs_wave_raw.value,
                     arr_obs_erro_raw.value,
                     np.zeros(len(arr_obs_flux_raw)), 
                     color='black', linewidth=0, alpha=.15, zorder=9)
    if templatename in ['vvds_elliptical','vvds_s0','red_galaxy']:
        alpha_abs = 1.0
        alpha_emi = 0.4
    else:
        alpha_abs = 0.6
        alpha_emi = 0.8
    xmin1 = arr_obs_wave_raw[ 0].value
    xmax1 = arr_obs_wave_raw[-1].value
    ymin1 = ax1.get_ylim()[0]
    ymax1 = ax1.get_ylim()[1]
    plot_lines(dic_emi_lines, dic_abs_lines, z_guess, ax1, 
               alpha_abs, alpha_emi, xmin1, xmax1, ymin1, ymax1)
    plot_blocks(ax1, xmin1, xmax1, ymin1, ymax1)
    plot_nega_flux(ax1, arr_obs_wave_raw, arr_obs_flux_raw)
    ax1.set_ylabel('Flux')
    ax1.minorticks_on()
    ax1.tick_params(  which='minor', length=3,bottom=True,top=True,direction='in')
    ax1.tick_params(  which='major', length=6,bottom=True,top=True,direction='in')
    ax1.set_xlabel(r'Observed Wavelength ($\AA$)')
    ax1.set_xlim(xmin1, xmax1)
    legend1 = ax1.legend(loc='upper left', prop={'size': 8})
    legend1.get_frame().set_facecolor('white')
    
    # Sub-plot for subtracted spectra
    ax4.plot(arr_obs_wave,
             arr_obs_rela_flux, 
             label='Observed', 
             color='black', linewidth=1, zorder=10)
    ax4.plot(arr_tem_wave     , 
             arr_tem_rela_flux, 
             label='Template', 
             color='red', linewidth=1, zorder=9, linestyle='-')
    xmin4 = xmin1
    xmax4 = xmax1
    ymin4 = np.min(arr_obs_rela_flux)
    ymax4 = np.min([1.5*np.max(arr_obs_rela_flux), 5*np.std(arr_obs_rela_flux)])
    plot_lines(dic_emi_lines, dic_abs_lines, z_guess, ax4, 
               alpha_abs, alpha_emi, xmin1, xmax1, ymin4, ymax4)
    plot_blocks(ax4, xmin4, xmax4, ymin4, ymax4)
    plot_nega_flux_rela(ax4, arr_obs_wave_raw, arr_obs_flux_raw, ymin4)
    ax4.set_ylabel('Relative Flux'+'\n'+' (spec / cont - 1)')
    ax4.minorticks_on()
    ax4.tick_params(which='minor', length=3,bottom=True,top=True,direction='in')
    ax4.tick_params(which='major', length=6,bottom=True,top=True,direction='in',
                    axis='x', labeltop=True, labelbottom=True)
    ax4.set_xlabel(r'Observed Wavelength ($\AA$)')
    ax4.set_xlim(xmin1, xmax1)
    ax4.set_ylim(ymin4, ymax4)
    legend4 = ax4.legend(loc='upper left', prop={'size': 8})
    legend4.get_frame().set_facecolor('white')
    
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
                       ls='-',  zorder=6, color='red', alpha=0.5, linewidth=4)
            if z_guess >= xmin2 and z_guess <= xmax2:
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
    except ValueError:
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
    return