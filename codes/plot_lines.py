def plot_lines(dic_emi_lines, dic_abs_lines, z_guess, ax0, 
               alpha_abs, alpha_emi, xmin1, xmax1, ymin1, ymax1):
    for key, values in dic_emi_lines.items():
        if key == '[OII]':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='blue', 
                           alpha=alpha_emi, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  "[O II]", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='blue', alpha=alpha_emi)
        if key == '[OIII]':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='blue', 
                           alpha=alpha_emi, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  "[O III]", size=10, zorder=9.1, 
                         horizontalalignment='center', verticalalignment='bottom', 
                         color='blue', alpha=alpha_emi)
        if key == '[SII]':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='blue', 
                           alpha=alpha_emi, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  "[S II]", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='blue', alpha=alpha_emi)
        if key == '[SIII]':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='blue', 
                           alpha=alpha_emi, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  "[S III]", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='blue', alpha=alpha_emi)
        if key == 'Ha':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='blue', 
                           alpha=alpha_emi, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  r"H$\alpha$", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='blue', alpha=alpha_emi)
        if key == 'Hb':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='blue', 
                           alpha=alpha_emi, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  r"H$\beta$", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='blue', alpha=alpha_emi)
        if key == 'Hg':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='blue', 
                           alpha=alpha_emi, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  r"H$\gamma$", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='blue', alpha=alpha_emi)
        if key == 'NII':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='blue', 
                           alpha=alpha_emi, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  "N II", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='blue', alpha=alpha_emi)
    for key, values in dic_abs_lines.items():
        if key == 'CaH&K':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='orange', 
                           alpha=alpha_abs, linewidth=1)
            if PlotOrNot:
                ax0.text(wave,  ymin1,  "Ca H&K", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='orange', alpha=alpha_abs)
        if key == 'G-band':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='orange', 
                           alpha=alpha_abs, linewidth=1)
            if PlotOrNot:
                ax0.text(wave-10,  ymin1,  "G-band", size=9, zorder=9.1, 
                         horizontalalignment='right', verticalalignment='bottom', 
                         color='orange', alpha=alpha_abs)
        if key == 'MgI':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='orange', 
                           alpha=alpha_abs, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  "Mg I", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='orange', alpha=alpha_abs)
        if key == 'Na':
            PlotOrNot = False
            for wave in values:
                wave *= (1+z_guess)
                if wave < xmin1 or wave > xmax1:
                    PlotOrNot = False
                    break
                PlotOrNot = True
                ax0.vlines(x=wave, ymin=ymin1, ymax=ymax1, 
                           ls='--',  zorder=8, color='orange', 
                           alpha=alpha_abs, linewidth=1)
            if PlotOrNot:
                ax0.text(wave+10,  ymin1,  "Na", size=10, zorder=9.1, 
                         horizontalalignment='left', verticalalignment='bottom', 
                         color='orange', alpha=alpha_abs)
    return
