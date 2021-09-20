# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 10:22:05 2021

@author: Jing
"""
import numpy as np
import pandas as pd
import matplotlib
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from .ions import ions_WEIGHT, ions_CHARGE

def plot(df, 
         unit='mg/L', 
         figname='Gaillardet diagram', 
         figformat='jpg'):
    """Plot the Gaillardet diagram.
    
    Parameters
    ----------
    df : class:`pandas.DataFrame`
        Geochemical data to draw Gibbs diagram.
    unit : class:`string`
        The unit used in df. Currently only mg/L is supported. 
    figname : class:`string`
        A path or file name when saving the figure.
    figformat : class:`string`
        The file format, e.g. 'png', 'pdf', 'svg'
        
        
    References
    ----------
    .. [1] Gaillardet, J. et al. 1999.
           Global silicate weathering and CO consumption rates deduced
           from the chemistry of large rivers
           Chemical Geology 159, 3-30.
           https://doi.org/10.1016/S0009-2541(99)00031-5
    """
    # Basic data check 
    # -------------------------------------------------------------------------
    # Determine if the required geochemical parameters are defined. 
    if not {'Ca', 'Mg', 'Na', 'HCO3'}.issubset(df.columns):
        raise RuntimeError("""
        Gibbs diagram uses geochemical parameters Ca, Mg, Na, and HCO3.
        Confirm that these parameters are provided.""")
        
    # Determine if the provided unit is allowed
    ALLOWED_UNITS = ['mg/L']
    if unit not in ALLOWED_UNITS:
        raise RuntimeError("""
        Currently only mg/L is supported.
        Convert the unit if needed.""")
        
    # Convert mg/L to meq/L
    gmol = np.array([ions_WEIGHT['Ca'], 
                     ions_WEIGHT['Mg'], 
                     ions_WEIGHT['Na'], 
                     ions_WEIGHT['HCO3']])

    eqmol = np.array([ions_CHARGE['Ca'], 
                      ions_CHARGE['Mg'], 
                      ions_CHARGE['Na'], 
                      ions_CHARGE['HCO3']])

    tmpdf = df[['Ca', 'Mg', 'Na', 'HCO3']]
    dat = tmpdf.values
    
    meqL = (dat / abs(gmol)) * abs(eqmol)
    
    # Do the plot
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(12, 10))
    
    # Plot the scatters
    ax1 = fig.add_subplot(221)
    ax1.loglog()
        
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            plt.scatter(meqL[i, 0] / meqL[i, 2], meqL[i, 3] / meqL[i, 2], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        label=TmpLabel, 
                        edgecolors='black')
            
        except(ValueError):
            pass
        
    # Creat the legend
    plt.legend(loc='lower right', markerscale=1, frameon=False, 
               labelspacing=0.25, handletextpad=0.25)
    
    # Show horizontal line at 100
    ax1.axhline(y=100, linestyle=':', linewidth=1, color='k')

    # Add a rectangle for Evaporites
    rect = mpatches.Rectangle([0.12, 0.15], 0.15, 0.20, 
                              fc="w", ec='k', alpha=0.2, hatch='///')
    ax1.add_patch(rect)
    ax1.text(0.3, 0.2, 'Evaporites', fontsize=12)
    
    # Add an ellipse for Silicates
    axins = inset_axes(ax1, width="100%", height="100%", loc=3)
    ellipse = mpatches.Ellipse([0.2, 0.35], 0.20, 0.13, angle=45, 
                               fc="w", ec='k', alpha=0.2, hatch='\\\\\\')
    axins.add_patch(ellipse)
    axins.axis('off')
    ax1.text(0.15, 4, 'Silicates', fontsize=12)
    
    # Add an ellipse for Carbonates
    axins = inset_axes(ax1, width="100%", height="100%", loc=2)
    ellipse = mpatches.Ellipse([0.85, 0.90], 0.25, 0.1, angle=30,
                               fc="w", ec='k', alpha=0.2, hatch='++')
    axins.add_patch(ellipse)
    axins.axis('off')
    ax1.text(3, 70, 'Carbonates', fontsize=12)

    # Set the labels and ticks
    ax1.set_xlabel('$Ca^{2+}/Na^+$', fontsize=12)
    ax1.set_ylabel('$HCO_3^{-}/Na^+$', fontsize=12)
  
    minorticks_off()
    tick_params(which='major', direction='in', length=4, width=1.25)
    tick_params(which='minor', direction='in', length=2.5, width=1.25)
    
    labels = ax1.get_xticklabels() + ax1.get_yticklabels()
    [label.set_fontsize(10) for label in labels]
    
    ax1.spines['top'].set_linewidth(1.25)
    ax1.spines['top'].set_color('k')
    ax1.spines['bottom'].set_linewidth(1.25)
    ax1.spines['bottom'].set_color('k')
    ax1.spines['left'].set_linewidth(1.25)
    ax1.spines['left'].set_color('k')
    ax1.spines['right'].set_linewidth(1.25)
    ax1.spines['right'].set_color('k')
    
    ax1.set_xlim(0.1, 100)
    ax1.set_ylim(0.1, 250)

    # -------------------------------------------------------------------------
    ax2 = fig.add_subplot(222)
    ax2.loglog()
    # Plot the scatters
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            plt.scatter(meqL[i, 0] / meqL[i, 2], meqL[i, 1] / meqL[i, 2], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        #label=TmpLabel, 
                        edgecolors='black')
            
        except(ValueError):
            pass
        
    # Show horizontal line at 10
    ax2.axhline(y=10, linestyle=':', linewidth=1, color='k')

    # Add a rectangle for Evaporites
    rect = mpatches.Rectangle([0.12, 0.015], 0.15, 0.020, 
                              fc="w", ec='k', alpha=0.2, hatch='///')
    ax2.add_patch(rect)
    ax2.text(0.3, 0.02, 'Evaporites', fontsize=12)
    
    # Add an ellipse for Silicates
    axins = inset_axes(ax2, width="100%", height="100%", loc=3)
    ellipse = mpatches.Ellipse([0.20, 0.40], 0.20, 0.13, angle=45, 
                               fc="w", ec='k', alpha=0.2, hatch='\\\\\\')
    axins.add_patch(ellipse)
    axins.axis('off')
    ax2.text(0.15, 0.6, 'Silicates', fontsize=12)
    
    # Add an ellipse for Carbonates
    axins = inset_axes(ax2, width="100%", height="100%", loc=2)
    ellipse = mpatches.Ellipse([0.87, 0.90], 0.20, 0.17, angle=45, 
                               fc="w", ec='k', alpha=0.2, hatch='++')
    axins.add_patch(ellipse)
    axins.axis('off')
    ax2.text(4, 7, 'Carbonates', fontsize=12)
    
    # Set the lables and ticks
    ax2.set_xlabel('$Ca^{2+}/Na^+$', fontsize=12)
    ax2.set_ylabel('$Mg^{2+}/Na^+$', fontsize=12)
    
    labels = ax2.get_xticklabels() + ax2.get_yticklabels()
    [label.set_fontsize(10) for label in labels]
    
    ax2.set_xlim(0.1, 100)
    ax2.set_ylim(0.01, 25)
    
    minorticks_on()
    tick_params(which='major', direction='in', length=4, width=1.25)
    tick_params(which='minor', direction='in', length=2.5, width=1.25)
    
    ax2.spines['top'].set_linewidth(1.25)
    ax2.spines['top'].set_color('k')
    ax2.spines['bottom'].set_linewidth(1.25)
    ax2.spines['bottom'].set_color('k')
    ax2.spines['left'].set_linewidth(1.25)
    ax2.spines['left'].set_color('k')
    ax2.spines['right'].set_linewidth(1.25)
    ax2.spines['right'].set_color('k')
   
    # Save the figure
    plt.savefig(figname + '.' + figformat, format=figformat, 
                bbox_inches='tight', dpi=300)
    
    return

if __name__ == '__main__':
    data = {'Sample' : ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample5'],
            'Label'  : ['C1', 'C2', 'C2', 'C3', 'C4', 'C4'],
            'Color'  : ['red', 'blue', 'blue', 'yellow', 'yellow', 'green'],
            'Marker' : ['o', 'o', 'o', 'o', 'o', 'o'],
            'Size'   : [30, 30, 30, 30, 30, 30],
            'Alpha'  : [0.6, 0.6, 0.6, 0.6, 0.6, 0.6],
            'pH'     : [7.78, 7.78, 7.85, 7.61, 7.45, 7.45],
            'Ca'     : [205.2, 214.5, 268.7, 215.8, 227.4, 221.8],
            'Mg'     : [63.77, 66.67, 58.9, 65.57, 69.86, 67.97],
            'Na'     : [21.36, 22.55, 25.76, 23.45, 32.63, 36.53],
            'K'      : [1.32, 2.14, 3.78, 2.64, 1.52, 4.24],
            'HCO3'   : [584.5, 584.5, 571.7, 557.1, 426.2, 484.1],
            'CO3'    : [0, 0, 0, 0, 0, 0],
            'Cl'     : [55.89, 56.09, 42.53, 65.27, 63.77, 63.28],
            'SO4'    : [308.4, 310.4, 521, 359.2, 448.1, 449.1],
            'NO3'    : [15.64, 14.78, 12.67, 16.2, 17.81, 14.51],
            'TDS'    : [1258.6, 1274.2, 1507, 1307, 1289.3, 1344.1],
            }
    df = pd.DataFrame(data)
    plot(df, unit='mg/L', figname='Gaillardet diagram', figformat='jpg')