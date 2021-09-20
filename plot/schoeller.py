# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 13:11:29 2021

@author: Jing
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .ions import ions_WEIGHT, ions_CHARGE

def plot(df, 
         unit='mg/L', 
         figname='Schoeller diagram', 
         figformat='jpg'):
    """Plot the HFE-D  diagram.
    
    Parameters
    ----------
    df : class:`pandas.DataFrame`
        Geochemical data to draw HFE-D diagram.
    unit : class:`string`
        The unit used in df. Currently only mg/L is supported. 
    figname : class:`string`
        A path or file name when saving the figure.
    figformat : class:`string`
        The file format, e.g. 'png', 'pdf', 'svg'
        
        
    References
    ----------
    .. [1] Güler, et al. 2002.
           Evaluation of graphical and multivariate statistical methods for classification of water chemistry data
           Hydrogeology Journal 10(4):455-474
           https://doi.org/10.1007/s10040-002-0196-6
    """
    # Determine if the required geochemical parameters are defined. 
    if not {'Ca', 'Mg', 'Na', 'K', 'Cl', 'SO4', 'HCO3'}.issubset(df.columns):
        raise RuntimeError("""
        Schoeller diagram uses geochemical parameters Ca, Mg, Na, K, Cl, SO4, and HCO3.
        Confirm that these parameters are provided.""")
        
    # Determine if the provided unit is allowed.
    ALLOWED_UNITS = ['mg/L']
    if unit not in ALLOWED_UNITS:
        raise RuntimeError("""
        Currently only mg/L is supported.
        Convert the unit if needed.""")
    
    # Convert mg/L to meq/L
    # -------------------------------------------------------------------------
    gmol = np.array([ions_WEIGHT['Ca'], 
                     ions_WEIGHT['Mg'], 
                     ions_WEIGHT['Na'], 
                     ions_WEIGHT['K'], 
                     ions_WEIGHT['Cl'], 
                     ions_WEIGHT['SO4'],
                     ions_WEIGHT['HCO3']])

    eqmol = np.array([ions_CHARGE['Ca'], 
                      ions_CHARGE['Mg'], 
                      ions_CHARGE['Na'], 
                      ions_CHARGE['K'],
                      ions_CHARGE['Cl'],
                      ions_CHARGE['SO4'],
                      ions_CHARGE['HCO3']])

    tmpdf = df[['Ca', 'Mg', 'Na', 'K', 'Cl', 'SO4', 'HCO3']]
    dat = tmpdf.values
    
    meqL = (dat / abs(gmol)) * abs(eqmol)
    
    # Do the plot
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)
    ax.semilogy()
    
    # Plot the lines
    # -------------------------------------------------------------------------
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            ax.plot([1, 2, 3, 4, 5, 6, 7], meqL[i, :], 
                    marker=df.at[i, 'Marker'],
                    color=df.at[i, 'Color'], 
                    alpha=df.at[i, 'Alpha'],
                    label=TmpLabel) 
        except(ValueError):
                pass
            
    # Background settings
    ax.set_xticks([1, 2, 3, 4, 5, 6, 7])
    ax.set_xticklabels(['$Ca^{2+}$', '$Mg^{2+}$', '$Na^+$', '$K^+$', 
                        '$Cl_-$', '$SO_4^-$', '$HCO_3^-$'])
    ax.set_ylabel('$meq/L$', fontsize=12)
    
    # Set the limits
    ax.set_xlim([1, 7])
    ax.set_ylim([np.min(meqL) * 0.1, np.max(meqL) * 1.1])
    
    # Plot the vertical lines
    for xtick in [1, 2, 3, 4, 5, 6, 7]:
        plt.axvline(xtick, linewidth=1, color='grey', linestyle='dashed')
            
    # Creat the legend
    ax.legend(loc='best', markerscale=1, frameon=False, 
              labelspacing=0.25, handletextpad=0.25)
    
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
    plot(df, unit='mg/L', figname='Schoeller diagram', figformat='jpg')
    
    
    