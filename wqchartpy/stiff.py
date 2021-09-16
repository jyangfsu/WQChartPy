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

from ions import ions_WEIGHT, ions_CHARGE

def plot_Stiff(df, 
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
    .. [1] Gibbs, Ronald J. 1970.
           Mechanisms Controlling World Water Chemistry.
           Science 170, 978–988.
           https://doi.org/10.1126/science.170.3962.1088
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
            plt.text(-cat_max, 2.9, 'Na+K', fontsize=12, ha= 'right')
            plt.text(-cat_max, 1.9, 'Ca', fontsize=12, ha= 'right')
            plt.text(-cat_max, 1.0, 'Mg', fontsize=12, ha= 'right')
            
            plt.text(an_max, 2.9,'Cl',fontsize=12, ha= 'left')
            plt.text(an_max, 1.9,'HCO'+'$_{3}$',fontsize=12,ha= 'left')
            plt.text(an_max, 1.0,'SO'+'$_{4}$',fontsize=12,ha= 'left')

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
            ax.set_xlabel('Stiff diagram (meq/L)', fontsize=12)
                
        except(ValueError):
                pass
    
        # Save the figure
        plt.savefig(figname + '_' + str(df.at[i, 'Sample']) + '.' + figformat, format=figformat, 
                    bbox_inches='tight', dpi=300)
        
    return
    

df = pd.read_csv('data.csv', nrows =5)   
plot_Stiff(df, unit='mg/L', figname='Stiff diagram', figformat='pdf')