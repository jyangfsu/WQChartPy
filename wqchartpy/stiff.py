# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 16:38:48 2021

@author: Jing
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from pylab import *

from .ions import ions_WEIGHT, ions_CHARGE

# Define the plotting function
def plot(df, 
         unit='mg/L', 
         figname='Stiff diagram', 
         figformat='jpg'):
    """Plot the Piper diagram.
    
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
    .. [1] Stiff, H.A. 1951.
           The Interpretation of Chemical Water Analysis by Means of Patterns
           Journal of Petroleum Technology 3(10): 15-3
           https://doi.org/10.2118/951376-G
    """
    # Basic data check 
    # -------------------------------------------------------------------------
    # Determine if the required geochemical parameters are defined. 
    if not {'Sample', 'Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4'}.issubset(df.columns):
        raise RuntimeError("""
        Stiff diagram uses geochemical parameters Ca, Mg, Na, K, HCO3, Cl, and SO4.
        Also, Sample is requied to save the Stiff diagram to disk for each sample.
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
                     ions_WEIGHT['K'], 
                     ions_WEIGHT['HCO3'],
                     ions_WEIGHT['Cl'], 
                     ions_WEIGHT['SO4']])

    eqmol = np.array([ions_CHARGE['Ca'], 
                      ions_CHARGE['Mg'], 
                      ions_CHARGE['Na'], 
                      ions_CHARGE['K'], 
                      ions_CHARGE['HCO3'],  
                      ions_CHARGE['Cl'], 
                      ions_CHARGE['SO4']])

    tmpdf = df[['Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4']]
    dat = tmpdf.values
    
    meqL = (dat / abs(gmol)) * abs(eqmol)
    cat_max = np.max(np.array(((meqL[:, 2] + meqL[:, 3]), meqL[:, 0], meqL[:, 1])))
    an_max = np.max(meqL[:, 4:])
    
    # Plot the Stiff diagrams for each sample
    # -------------------------------------------------------------------------
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            #xtickpositions = [-10, -5, 0, 5, 10]
            
            x = [-(meqL[i, 2] + meqL[i, 3]), -meqL[i, 0], -meqL[i, 1], 
                 meqL[i, 6], meqL[i, 4], meqL[i, 5], -(meqL[i, 2] + meqL[i, 3])]
            y = [3, 2, 1, 1, 2, 3, 3]
            
            plt.figure(figsize=(3, 3))
            plt.fill(x, y, facecolor='w', edgecolor='k', linewidth=1.25)
            plt.plot([0, 0], [1, 3], 'k-.', linewidth=1.25)
            plt.plot([-0.5, 0.5], [2,2], 'k')
            plt.plot([-0.5, 0.5], [3,3], 'k')
            plt.xlim([-cat_max, an_max])
            plt.text(-cat_max, 2.9, 'Na$^+$' + '+' + 'K$^+$', fontsize=12, ha= 'right')
            plt.text(-cat_max, 1.9, 'Ca$^{2+}$', fontsize=12, ha= 'right')
            plt.text(-cat_max, 1.0, 'Mg$^{2+}$', fontsize=12, ha= 'right')
            
            plt.text(an_max, 2.9,'Cl$^-$',fontsize=12, ha= 'left')
            plt.text(an_max, 1.9,'HCO'+'$_{3}^-$',fontsize=12,ha= 'left')
            plt.text(an_max, 1.0,'SO'+'$_{4}^{2-}$',fontsize=12,ha= 'left')

            #plt.xticks(xtickpositions, ['10', '5', '0', '5', '10', '20']) 
            ax = plt.gca()
            ax.spines['left'].set_color('None')
            ax.spines['right'].set_color('None')
            ax.spines['top'].set_color('None')
            minorticks_off()
            tick_params(which='major', direction='out', length=4, width=1.25)
            tick_params(which='minor', direction='in', length=2, width=1.25)
            ax.spines['bottom'].set_linewidth(1.25)
            ax.spines['bottom'].set_color('k')
            #ylim(0.8, 3.2)
            setp(gca(), yticks=[], yticklabels=[])
            labels = ax.get_xticklabels()
            [label.set_fontsize(10) for label in labels]
            ax.set_xlabel('Stiff diagram (meq/L)', fontsize=12, weight='normal')
                
        except(ValueError):
                pass
    
        # Display the info
        print("Stiff plot created for %s. Saving it now...\n" %str(df.at[i, 'Sample']))
    
        # Save the figure
        plt.savefig(figname + '_' + str(df.at[i, 'Sample']) + '.' + figformat, format=figformat, 
                    bbox_inches='tight', dpi=300)
        
    return

if __name__ == '__main__':
    # Example data
    data = {'Sample' : ['sample1', 'sample2', 'sample3', 'sample4', 'sample5', 'sample6'],
            'Label'  : ['C1', 'C2', 'C2', 'C3', 'C3', 'C1'],
            'Color'  : ['red', 'green', 'green', 'blue', 'blue', 'red'],
            'Marker' : ['o', 'o', 'o', 'o', 'o', 'o'],
            'Size'   : [30, 30, 30, 30, 30, 30],
            'Alpha'  : [0.6, 0.6, 0.6, 0.6, 0.6, 0.6],
            'pH'     : [7.8, 7.6, 7.5, 7.7, 7.4, 7.1],
            'Ca'     : [32, 46, 54, 50, 50, 134],
            'Mg'     : [6, 11, 11, 11, 22, 21],
            'Na'     : [28, 17, 16, 25, 25, 39],
            'K'      : [2.8, 0.7, 2.4, 2.8, 0.5, 6.4],
            'HCO3'   : [73, 201, 207, 244, 305, 275],
            'CO3'    : [0, 0, 0, 0, 0, 0],
            'Cl'     : [43, 14, 18, 18, 11, 96],
            'SO4'    : [48, 9, 10, 9, 9, 100],
            'TDS'    : [233, 299, 377, 360, 424, 673],
            }
    df = pd.DataFrame(data)
    # df = pd.read_csv('../data/data_template.csv')
    plot(df, unit='mg/L', figname='Stiff diagram', figformat='jpg')
    
    
    
    