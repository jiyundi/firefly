import numpy as np
from codes.plot_lines     import plot_lines
from codes.plot_blocks    import plot_blocks, plot_nega_flux, plot_nega_flux_rela

def plots_spec_z_combo(ax_2dspecimg,
                       farr, vmin, vmax, warr, spec2dbox_ratio,
                       ax_orig_spec, ax_rela_spec, 
                       ax_corr_dist, ax_peak_dist,
                       ax_orig_spe2, ax_rela_spe2,
                       dic_abs_lines, dic_emi_lines,
                       spec1dobjname, templatename, z_peak, z_guess,
                       lag, corr_normalized,
                       arr_obs,           # Figure 1+3: 1/4=arr_obs
                       arr_obs_cont,      # Figure 1+3: 2/4 (1/4=arr_obs)
                       arr_tem1_masked,   # Figure 1+3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
                       arr_tem1cont,      # Figure 1+3, 4/4
                       arr_obs_rela,    # Figure 2+4, 1/2
                       arr_tem1rela,    # Figure 2+4, 2/2
                       arr_tem1_masked_newz,     # Figure 3: 3/4 (Can't use arr_tem1 b/c un-normalized!)
                       arr_tem1cont_masked_newz, # Figure 3, 4/4
                       arr_tem1relaC,            # Figure 4, 2/2
                       colors):
    
    
    if templatename in ['vvds_elliptical','vvds_s0','red_galaxy',
                        'sdss_early_type','sdss_luminous_red']:
        alpha_abs = 1.0
        alpha_emi = 0.4
    else:
        alpha_abs = 0.6
        alpha_emi = 0.8
    
    # -------------------------------------------------------------------------------------
    # Sub-plot for observed and template spectra
    ax_orig_spec.plot(arr_obs[0], # Figure 1+3: 1/4
                      arr_obs[1], 
                      label='Observed spectrum',
                      color='black', linewidth=1, zorder=11)
    ax_orig_spec.plot(arr_obs_cont[0], # Figure 1+3: 2/4
                      arr_obs_cont[1], 
                      label='Observed continuum', alpha=0.75,
                      color='green', linewidth=2, zorder=10, linestyle=':')
    ax_orig_spec.plot(arr_tem1_masked[0], # Figure 1+3: 3/4
                      arr_tem1_masked[1], 
                      label='Template at z='+"{:.4f}".format(z_peak),
                      color=colors[0], linewidth=1,  zorder=11, linestyle='-')
    ax_orig_spec.plot(arr_tem1cont[0], # Figure 1+3: 4/4
                      arr_tem1cont[1], 
                      label='Template continuum', alpha=0.75, 
                      color=colors[0], linewidth=2, zorder=10, linestyle=':')
    ax_orig_spec.fill_between(arr_obs[0], # Figure 1+3: 1/4 (erro)
                              arr_obs[2],
                              np.zeros(len(arr_obs[1])), 
                              color='black', linewidth=0, alpha=.15, zorder=9)
    xmin1          = arr_obs[0][ 0]
    xmax_orig_spec = arr_obs[0][-1]
    ymin1          = ax_orig_spec.get_ylim()[0]
    ymax_orig_spec = ax_orig_spec.get_ylim()[1]
    plot_lines(dic_emi_lines, dic_abs_lines, z_peak, ax_orig_spec, 
               alpha_abs, alpha_emi, xmin1, xmax_orig_spec, ymin1, ymax_orig_spec)
    plot_blocks(ax_orig_spec, xmin1, xmax_orig_spec, ymin1, ymax_orig_spec)
    plot_nega_flux(ax_orig_spec, arr_obs[0], arr_obs[1])
    ax_orig_spec.set_ylabel('Flux')
    ax_orig_spec.minorticks_on()
    ax_orig_spec.tick_params(  which='minor', length=6,bottom=True,top=True,direction='in')
    ax_orig_spec.tick_params(  which='major', length=9,bottom=True,top=True,direction='in')
    # ax_orig_spec.set_xlabel(r'Observed Wavelength ($\AA$)', labelpad=-2)
    ax_orig_spec.set_xlim(xmin1, xmax_orig_spec)
    ax_orig_spec.set_ylim(ymin1, np.max(arr_obs[1]))
    legend1 = ax_orig_spec.legend(loc='upper left', prop={'size': 8})
    legend1.get_frame().set_facecolor('white')
    
    
    # -------------------------------------------------------------------------------------
    # plot 2d spectrum (mask reference row: 10)
    i_xmin1          = np.argmax(warr[10][warr[10] <= xmin1         ]) 
    i_xmax_orig_spec = np.argmax(warr[10][warr[10] <= xmax_orig_spec])
    farr_draw = farr[:, i_xmin1:i_xmax_orig_spec]
    ax_2dspecimg.imshow(farr_draw, 
                        cmap='gray', 
                        vmin=vmin, vmax=vmax,
                        extent=[xmin1, xmax_orig_spec, 0, len(farr[:,0])],  
                        # extent=[left, right, bottom, top]
                        aspect=np.sqrt((len(farr_draw[0])/len(farr_draw)) / spec2dbox_ratio)
                        )
    
    
    # -------------------------------------------------------------------------------------
    # Sub-plot for subtracted spectra
    ax_rela_spec.plot(arr_obs_rela[0], # Figure 2+4: 1/2
                      arr_obs_rela[1], 
                      label='Observed', 
                      color='black', linewidth=1, zorder=10)
    ax_rela_spec.plot(arr_tem1rela[0], # Figure 2+4: 2/2
                      arr_tem1rela[1], 
                      label='Template', 
                      color=colors[0], linewidth=1, zorder=9, linestyle='-')
    xmin4 = xmin1
    xmax_rela_spec = xmax_orig_spec
    ymin4 = np.min(arr_obs_rela[1])
    ymax_rela_spec = np.min([1.5*np.max(arr_obs_rela[1]), 
                             8*np.std(arr_obs_rela[1])])
    plot_lines(dic_emi_lines, dic_abs_lines, z_peak, ax_rela_spec, 
               alpha_abs, alpha_emi, xmin1, xmax_orig_spec, ymin4, ymax_rela_spec)
    plot_blocks(ax_rela_spec, xmin4, xmax_rela_spec, ymin4, ymax_rela_spec)
    plot_nega_flux_rela(ax_rela_spec, arr_obs[0], arr_obs[1], ymin4)
    ax_rela_spec.set_ylabel('Relative Flux'+'\n'+' (spec / cont - 1)')
    ax_rela_spec.minorticks_on()
    ax_rela_spec.tick_params(which='minor', length=3,bottom=True,top=True,direction='in')
    ax_rela_spec.tick_params(which='major', length=6,bottom=True,top=True,direction='in',
                    axis='x', labeltop=False, labelbottom=True)
    # ax_rela_spec.set_xlabel(r'Observed Wavelength ($\AA$)', labelpad=-2)
    ax_rela_spec.set_xlim(xmin1, xmax_orig_spec)
    ax_rela_spec.set_ylim(ymin4, ymax_rela_spec)
    legend4 = ax_rela_spec.legend(loc='upper left', prop={'size': 8})
    legend4.get_frame().set_facecolor('white')
    
    ax_orig_spec.text(0.5, 0.5, "{:.4f}".format(z_peak), size=60,
             horizontalalignment='center', verticalalignment='center', 
             color='olive', alpha=0.15, zorder=9,
             transform=ax_orig_spec.transAxes)
    ax_rela_spec.text(0.5, 0.5, "{:.4f}".format(z_peak), size=60,
             horizontalalignment='center', verticalalignment='center', 
             color='olive', alpha=0.15, zorder=9,
             transform=ax_rela_spec.transAxes)
    
    # -------------------------------------------------------------------------------------
    ax_orig_spe2.plot(arr_obs[0], # Figure 1+3: 1/4
                      arr_obs[1], 
                      label='Observed spectrum',
                      color='black', linewidth=1, zorder=11)
    ax_orig_spe2.plot(arr_obs_cont[0], # Figure 1+3: 2/4
                      arr_obs_cont[1], 
                      label='Observed continuum', alpha=0.75,
                      color='green', linewidth=2, zorder=10, linestyle=':')
    ax_orig_spe2.fill_between(arr_obs[0], # Figure 1+3: 1/4 (erro)
                              arr_obs[2],
                              np.zeros(len(arr_obs[1])), 
                              color='black', linewidth=0, alpha=.15, zorder=9)
    xmin1          = arr_obs[0][ 0]
    xmax_orig_spe2 = arr_obs[0][-1]
    ymin1          = ax_orig_spe2.get_ylim()[0]
    ymax_orig_spe2 = ax_orig_spe2.get_ylim()[1]
    plot_blocks(ax_orig_spe2, xmin1, xmax_orig_spe2, ymin1, ymax_orig_spe2)
    plot_nega_flux(ax_orig_spe2, arr_obs[0], arr_obs[1])
    ax_orig_spe2.set_ylabel('Flux')
    ax_orig_spe2.minorticks_on()
    ax_orig_spe2.tick_params(  which='minor', length=6,bottom=True,top=True,direction='in')
    ax_orig_spe2.tick_params(  which='major', length=9,bottom=True,top=True,direction='in')
    # ax_orig_spe2.set_xlabel(r'Observed Wavelength ($\AA$)', labelpad=-2)
    ax_orig_spe2.set_xlim(xmin1, xmax_orig_spe2)
    ax_orig_spe2.set_ylim(ymin1, np.max(arr_obs[1]))
    
    ax_rela_spe2.plot(arr_obs_rela[0], # Figure 2+4: 1/2
                      arr_obs_rela[1],
                      label='Observed', 
                      color='black', linewidth=1, zorder=10)
    xmin4 = xmin1
    xmax_rela_spe2 = xmax_orig_spe2
    ymin4 = np.min(arr_obs_rela[1])
    ymax_rela_spe2 = np.min([1.5*np.max(arr_obs_rela[1]), 
                             8*np.std(arr_obs_rela[1])])
    plot_blocks(ax_rela_spe2, xmin4, xmax_rela_spe2, ymin4, ymax_rela_spe2)
    plot_nega_flux_rela(ax_rela_spe2, arr_obs[0], arr_obs[1], ymin4)

    if z_guess != z_peak:
        ax_orig_spe2.plot(arr_tem1_masked_newz[0], # Figure 3: 3/4
                          arr_tem1_masked_newz[1], 
                          label='Template at z='+"{:.4f}".format(z_guess),
                          color=colors[1], linewidth=1,  
                          zorder=11, linestyle='-')
        ax_orig_spe2.plot(arr_tem1cont_masked_newz[0], # Figure 3: 4/4
                          arr_tem1cont_masked_newz[1],
                          label='Template continuum',
                          color=colors[1], linewidth=2, 
                          zorder=9, linestyle=':')
        ax_rela_spe2.plot(arr_tem1relaC[0], # Figure 4: 2/2
                          arr_tem1relaC[1], 
                          label='Template', alpha=0.75,
                          color=colors[1], linewidth=1, 
                          zorder=10, linestyle='-')
        
        ax_orig_spe2.text(0.5, 0.5, "{:.4f}".format(z_guess), size=60,
                          horizontalalignment='center', verticalalignment='center', 
                          color='purple', alpha=0.15, zorder=9,
                          transform=ax_orig_spe2.transAxes)
        ax_rela_spe2.text(0.5, 0.5, "{:.4f}".format(z_guess), size=60,
                          horizontalalignment='center', verticalalignment='center', 
                          color='purple', alpha=0.15, zorder=9,
                          transform=ax_rela_spe2.transAxes)
        
    plot_lines(dic_emi_lines, dic_abs_lines, z_guess, ax_orig_spe2, 
                alpha_abs, alpha_emi, xmin1, xmax_orig_spe2, ymin1, ymax_orig_spe2)
    plot_lines(dic_emi_lines, dic_abs_lines, z_guess, ax_rela_spe2, 
                alpha_abs, alpha_emi, xmin1, xmax_orig_spe2, ymin4, ymax_rela_spe2)
    
    ax_rela_spe2.set_ylabel('Relative Flux'+'\n'+' (spec / cont - 1)')
    ax_rela_spe2.minorticks_on()
    ax_rela_spe2.tick_params(which='minor', length=3,bottom=True,top=True,direction='in')
    ax_rela_spe2.tick_params(which='major', length=6,bottom=True,top=True,direction='in',
                    axis='x', labeltop=False, labelbottom=True)
    ax_rela_spe2.set_xlim(xmin1, xmax_orig_spe2)
    ax_rela_spe2.set_ylim(ymin4, ymax_rela_spe2)
    legend5 = ax_orig_spe2.legend(loc='upper left', prop={'size': 8})
    legend6 = ax_rela_spe2.legend(loc='upper left', prop={'size': 8})
    legend5.get_frame().set_facecolor('white')
    legend6.get_frame().set_facecolor('white')
    
    
    # -------------------------------------------------------------------------------------
    # Sub-plot for redshift probability distribution and guess
    ax_corr_dist.grid(alpha=0.75, linestyle='--', zorder=-1)
    ax_corr_dist.plot(lag, 
            corr_normalized, 
            color='green', linewidth=1, zorder=10)
    ax_corr_dist.set_xlim(left=0, right=None)
    ax_corr_dist.set_ylim(bottom=0, top=None)
    
    
    # Sub-plot for zoom-in redshift peak
    try:
        ax_peak_dist.grid(alpha=0.5, linestyle='--', zorder=-1)
        ax_peak_dist.plot(lag, 
                corr_normalized, 
                color='green', linewidth=1, zorder=10)
        z_peak = lag[np.argmax(corr_normalized)]
        ymin  = np.min(corr_normalized)
        ymax  = np.max(corr_normalized)*1.25
        xmin2 = np.min(lag[corr_normalized > np.max(corr_normalized)*0.5]) - 0.02
        ymin2 = np.max(corr_normalized)*0.0
        xmax_corr_peak_dist = np.max(lag[corr_normalized > np.max(corr_normalized)*0.5]) + 0.02
        ymax_corr_peak_dist = np.max(corr_normalized)*1.25
        ax_corr_dist.vlines(x=z_peak,  
                            ymin=ymin, 
                            ymax=ymax, 
                            ls='-', zorder=8, 
                            color='olive', alpha=0.5, linewidth=8)
        ax_peak_dist.text(z_peak, ymax_corr_peak_dist, 
                          r"$z_\mathrm{peak}=$"+"{:.4f}".format(z_peak), size=10,
                          horizontalalignment='center', verticalalignment='top', 
                          color='olive', zorder=7)
        ax_peak_dist.vlines(x=z_peak,
                            ymin=ymin2,
                            ymax=ymin2+(ymax_corr_peak_dist-ymin2)*0.9,
                            ls='-',  zorder=6,
                            color='olive', alpha=0.5, linewidth=4)
            
        if z_guess != z_peak:
            ax_corr_dist.vlines(x=z_guess, 
                                ymin=ymin, 
                                ymax=ymax, 
                                ls='-',  zorder=6, 
                                color='purple', alpha=0.3, linewidth=8)
            xmin2               = np.min([z_peak-0.05, z_guess-0.05])
            xmax_corr_peak_dist = np.max([z_peak+0.05, z_guess+0.05])
            ax_peak_dist.vlines(x=z_guess, 
                                ymin=ymin2, 
                                ymax=ymin2+(ymax_corr_peak_dist-ymin2)*0.8, 
                                ls='-',  zorder=8, 
                                color='purple', alpha=0.3, linewidth=4)
            ax_peak_dist.text(z_guess,  ymin2+(ymax_corr_peak_dist-ymin2)*0.9,  
                              r"$z_\mathrm{guess}=$"+"{:.4f}".format(z_guess), size=10, 
                              horizontalalignment='center', verticalalignment='top', 
                              color='purple', zorder=9)
        ax_peak_dist.set_xlim(xmin2, xmax_corr_peak_dist)
        ax_peak_dist.set_ylim(ymin2, ymax_corr_peak_dist)
    except ValueError:
        ax_peak_dist.text(0.5, 0.5, 'Zoom-in Unavailable', size=20,
                 horizontalalignment='center', verticalalignment='center', 
                 color='red', zorder=7,
                 transform=ax_peak_dist.transAxes)
        print('[Error] Zoom-in failed: check if any index is invalid!')
        pass
    ax_corr_dist.set_ylabel('Correlation')
    ax_peak_dist.set_ylabel('Correlation')
    ax_corr_dist.minorticks_on()
    ax_peak_dist.minorticks_on()
    ax_corr_dist.tick_params(which='minor', length=3, bottom=True,top=True, direction='in')
    ax_corr_dist.tick_params(which='major', length=9, bottom=True,top=True, direction='in',
                    axis='x', labeltop=False, labelbottom='on')
    ax_peak_dist.tick_params(which='minor', length=3, bottom=True,top=True, direction='in')
    ax_peak_dist.tick_params(which='major', length=15,bottom=True,top=True, direction='in',
                    axis='x', labeltop=False, labelbottom='on', labelsize=10)
    ax_corr_dist.set_xlabel(r'Redshift ($z$)')
    ax_peak_dist.set_xlabel(r'Redshift ($z$)')
    return