import numpy as np
from matplotlib.patches import Rectangle

def plot_blocks(ax0, xmin, xmax, ymin, ymax):
    # A-band (atm absorption)
    aband   = ax0.add_patch(Rectangle((7586, ymin), 
                            7708-7586, ymax-ymin, 
                            facecolor=(0.1,0.3,0.6,0.1), 
                            edgecolor=(0.0,0.0,0.0,0.2),
                            linewidth=1, zorder=1))
    aband.set_hatch('\\\\')
    
    # B-band (atm absorption)
    bband   = ax0.add_patch(Rectangle((6864, ymin), 
                            6945-6864, ymax-ymin, 
                            facecolor=(0.1,0.3,0.6,0.1), 
                            edgecolor=(0.0,0.0,0.0,0.2),
                            linewidth=1, zorder=1))
    bband.set_hatch('\\\\')
    
    # Edges
    edgeL   = ax0.add_patch(Rectangle((xmin, ymin), 
                            0.03*(xmax-xmin), ymax-ymin, 
                            facecolor=(0.1,0.1,0.1,0.2), 
                            edgecolor=(0.0,0.0,0.0,0.4),
                            linewidth=1, zorder=1))
    edgeR   = ax0.add_patch(Rectangle((xmax-0.03*(xmax-xmin), ymin), 
                            0.03*(xmax-xmin), ymax-ymin, 
                            facecolor=(0.1,0.1,0.1,0.2), 
                            edgecolor=(0.0,0.0,0.0,0.4),
                            linewidth=1, zorder=1))
    edgeL.set_hatch('///')
    edgeR.set_hatch('///')
    return

def plot_nega_flux(ax0, arr_obs_wave_raw, arr_obs_flux_raw):
    # Negative flux warning
    ymin = np.min(arr_obs_flux_raw)
    if ymin <= 0:
        mask = arr_obs_flux_raw > 0
        split_start = []
        split_end   = []
        n = 0
        for i in range(len(mask)-1):
            if mask[i]==False and n==0:
                split_start.append(i)
                if mask[i+1]==True:
                    split_end.append(i+1)
                    n = 0
                if mask[i+1]==False:
                    n += 1
            elif mask[i]==False and n!=0:
                if mask[i+1]==True:
                    split_end.append(i+1)
                    n = 0
                if mask[i+1]==False:
                    n += 1
        for i in range(len(split_start)):
            try:
                x1   = arr_obs_wave_raw[split_start[i]].value
                base = arr_obs_wave_raw[split_end[i]  ].value - x1
                wrong_f = ax0.add_patch(Rectangle((
                                        x1, ymin), 
                                        base, 0-ymin, 
                                        facecolor=(1.0,0.0,0.0,0.1), 
                                        edgecolor=(0.5,0.0,0.0,0.3),
                                        linewidth=1, zorder=1))
                wrong_f.set_hatch('xxxx')
            except IndexError:
                print('IndexError: list index out of range (occured in function "plot_nega_flux")')
    return

def plot_nega_flux_rela(ax0, arr_obs_wave_raw, arr_obs_flux_raw, ymin4):
    # Negative flux warning
    ymin = np.min(arr_obs_flux_raw)
    ymax = np.max(arr_obs_flux_raw)
    if ymin <= 0:
        mask = arr_obs_flux_raw > 0
        split_start = []
        split_end   = []
        n = 0
        for i in range(len(mask)-1):
            if mask[i]==False and n==0:
                split_start.append(i)
                if mask[i+1]==True:
                    split_end.append(i+1)
                    n = 0
                if mask[i+1]==False:
                    n += 1
            elif mask[i]==False and n!=0:
                if mask[i+1]==True:
                    split_end.append(i+1)
                    n = 0
                if mask[i+1]==False:
                    n += 1
        for i in range(len(split_start)):
            try:
                x1   = arr_obs_wave_raw[split_start[i]].value
                base = arr_obs_wave_raw[split_end[i]  ].value - x1
                wrong_f = ax0.add_patch(Rectangle((
                                        x1, ymin), 
                                        base, ymax-ymin, 
                                        facecolor=(1.0,0.0,0.0,0.1), 
                                        edgecolor=(0.5,0.0,0.0,0.3),
                                        linewidth=1, zorder=1))
                wrong_f.set_hatch('xxxx')
            except IndexError:
                print('IndexError: list index out of range (occured in function "plot_nega_flux_rela")')
    return
