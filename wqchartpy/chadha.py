# -*- coding: utf-8 -*-
"""
Created on Thu Sep 16 11:36:50 2021

@author: Jing
"""
# Import modules
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .ions import ions_WEIGHT, ions_CHARGE

# Define the Chadha plotting function
def plot(df, 
         unit='mg/L', 
         figname='Chadha diagram', 
         figformat='jpg'):
    """Plot the Chadha diagram.
    
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
    .. [1] Chadha, D. K. 1999.
           A proposed new diagram for geochemical classification of natural waters and interpretation of chemical data
           Hydrogeology Journal 7:431–439
           https://doi.org/10.1007/s100400050216
    """
    # Determine if the required geochemical parameters are defined. 
    if not {'Ca', 'Mg', 'Na', 'K', 'HCO3', 'CO3', 'Cl', 'SO4'}.issubset(df.columns):
        raise RuntimeError("""
        Chada diagram uses geochemical parameters Ca, Mg, Na, K, HCO3, CO3, Cl, and SO4.
        Confirm that these parameters are provided.""")
        
    # Determine if the provided unit is allowed.
    ALLOWED_UNITS = ['mg/L']
    if unit not in ALLOWED_UNITS:
        raise RuntimeError("""
        Currently only mg/L is supported.
        Convert the unit if needed.""")
        
    # Change default settings for figures
    # -------------------------------------------------------------------------
    plt.style.use('default')
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.labelsize'] = 10
    plt.rcParams['axes.labelweight'] = 'bold'
    plt.rcParams['axes.titlesize'] = 10
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.rcParams['legend.fontsize'] = 10
    plt.rcParams['figure.titlesize'] = 10
    
    # Plot background
    # -------------------------------------------------------------------------
    xmin = -100
    xmax = +100
    ymin = -100
    ymax = +100
    
    plt.figure(figsize=(8, 8))
    ax = plt.subplot(111)
    
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.spines['left'].set_position(('data', 0))
    ax.spines['bottom'].set_position(('data', 0))
    
    # Plot the domain
    plt.plot([xmin, xmax, xmax, xmin, xmin], [ymin, ymin, ymax, ymax, ymin],
             linestyle='-', linewidth=1.5, color='k')
    
    plt.xticks([-80, -60, -40, -20, 20, 40, 60, 80],
               ['-80', '-60', '-40', '-20', '+20', '+40', '+60', '+80'])
    plt.yticks([-80, -60, -40, -20, 20, 40, 60, 80],
               ['-80', '-60', '-40', '-20', '+20', '+40', '+60', '+80'])
    
    # The axis lines in the center
    ax.plot([xmin, xmax], [0, 0], linestyle='-', linewidth=1.0, color='k')
    ax.plot([0, 0], [ymin, ymax], linestyle='-', linewidth=1.0, color='k')
    
    # Labels
    plt.text(0, ymin * 1.08, '$(Ca^{2+}+Mg^{2+})-(Na^++K^+)$' + '\nMilliequivalent percentage', 
             ha='center', va='center', fontsize=12)
    plt.arrow(xmin * 0.7, ymin * 1.05, 0.1 * (xmax- xmin), 0, fc='k', head_width=2, head_length=4)
    plt.arrow(47, -105, 20, 0, fc='k', head_width=2, head_length=4)
    
    plt.text(-107.5, 0, 'Milliequivalent percentage\n' + '$(HCO_3^-+CO_3^{2-})-(Cl^-+SO_4^{2-})$', 
             ha='center', va='center', fontsize=12, rotation=90)
    plt.arrow(-105, -75, 0, 20, fc='k', head_width=2, head_length=4)
    plt.arrow(-105, 50, 0, 20, fc='k', head_width=2, head_length=4)
    
    # Mark the eight subfileds
    plt.text(50, 0, '1', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(-50, 0, '2', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(0, 50, '3', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(0, -50, '4', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(50, 50, '5', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(50, -50, '6', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(-50, -50, '7', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(-50, 50, '8', fontsize=26, color="0.6", 
             ha='center', va='center')
    
    
    
    
    # Convert mg/L to meq/L
    # -------------------------------------------------------------------------
    gmol = np.array([ions_WEIGHT['Ca'], 
                     ions_WEIGHT['Mg'], 
                     ions_WEIGHT['Na'], 
                     ions_WEIGHT['K'], 
                     ions_WEIGHT['HCO3'],
                     ions_WEIGHT['CO3'], 
                     ions_WEIGHT['Cl'], 
                     ions_WEIGHT['SO4']])

    eqmol = np.array([ions_CHARGE['Ca'], 
                      ions_CHARGE['Mg'], 
                      ions_CHARGE['Na'], 
                      ions_CHARGE['K'], 
                      ions_CHARGE['HCO3'], 
                      ions_CHARGE['CO3'], 
                      ions_CHARGE['Cl'], 
                      ions_CHARGE['SO4']])

    tmpdf = df[['Ca', 'Mg', 'Na', 'K', 'HCO3', 'CO3', 'Cl', 'SO4']]
    dat = tmpdf.values
    
    meqL = (dat / abs(gmol)) * abs(eqmol)
    
    # Calculate the percentages
    # -------------------------------------------------------------------------
    sumcat = np.sum(meqL[:, 0:4], axis=1)
    suman = np.sum(meqL[:, 4:], axis=1)
    cat = np.zeros((dat.shape[0], 3))
    an = np.zeros((dat.shape[0], 3))
    cat[:, 0] = meqL[:, 0] / sumcat                  # Ca
    cat[:, 1] = meqL[:, 1] / sumcat                  # Mg
    cat[:, 2] = (meqL[:, 2] + meqL[:, 3]) / sumcat   # Na+K
    an[:, 0] = (meqL[:, 4] + meqL[:, 5]) / suman     # HCO3 + CO3
    an[:, 2] = meqL[:, 6] / suman                    # Cl
    an[:, 1] = meqL[:, 7] / suman                    # SO4
    
    # Plot the scatter
    # -------------------------------------------------------------------------
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            ax.scatter(100 * (cat[i, 0] +  cat[i, 1] - cat[i, 2]), 
                       100 * (an[i, 0] - (an[i, 1] + an[i, 2])), 
                       marker=df.at[i, 'Marker'],
                       s=df.at[i, 'Size'], 
                       color=df.at[i, 'Color'], 
                       alpha=df.at[i, 'Alpha'],
                       label=TmpLabel, 
                       edgecolors='black') 
        except(ValueError):
            pass
            
    # Creat the legend
    ax.legend(bbox_to_anchor=(0.085, 0.95), markerscale=1, fontsize=12,
              frameon=False, 
              labelspacing=0.25, handletextpad=0.25)
    
    # Display the info
    print("Chada plot created. Saving it now...\n")
    
    # Save the figure
    
    plt.savefig(figname + '.' + figformat, format=figformat, 
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
    plot(df, unit='mg/L', figname='Chada diagram', figformat='jpg')

    
