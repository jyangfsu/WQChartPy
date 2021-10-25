# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 10:36:15 2020

@author: Jing
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .ions import ions_WEIGHT, ions_CHARGE

# Define the plotting function
def plot(df, 
         unit='mg/L', 
         figname='HFE-D diagram', 
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
    .. [1] Gim´enez-Forcada, E. 2009.
           Dynamic of Sea Water Interface using Hydrochemical Facies Evolution Diagram
           Groundawter 48(2), 212–216.
           https://doi.org/10.1111/j.1745-6584.2009.00649.x
    """
    # Determine if the required geochemical parameters are defined. 
    if not {'Ca', 'Mg', 'Na', 'K', 'HCO3', 'CO3', 'Cl', 'SO4'}.issubset(df.columns):
        raise RuntimeError("""
        HFE-D uses geochemical parameters Ca, Mg, Na, K, HCO3, CO3, Cl, and SO4.
        Confirm that these parameters are provided.""")
        
    # Determine if the provided unit is allowed.
    ALLOWED_UNITS = ['mg/L']
    if unit not in ALLOWED_UNITS:
        raise RuntimeError("""
        Currently only mg/L is supported.
        Convert the unit if needed.""")
    '''
    # Seawater concentrations from Turekian, K.K. ,1968.- Oceans , Prentice Hall
    CONC_SWAWATER = {'Ca': 411,
                     'Mg': 1290,
                     'Na': 10800,
                     'K': 392,
                     'HCO3': 142,
                     'SO4': 2712,
                     'Cl': 19400}
    '''
    # Seawater concentrations from Table 2.3 of Appelo and Postma, 2005.- Geochemistry, Groundwater and Pollution
    CONC_SWAWATER = {'Ca': 429,      # 10.7  mmol/L
                     'Mg': 1339,     # 55.1  mmol/L
                     'Na': 11150,    # 485.0 mmol/L
                     'K': 414,       # 10.6  mmol/L
                     'HCO3': 146,    # 2.4   mmol/L 
                     'SO4': 2815,    # 29.3  mmol/L
                     'Cl': 20066}    # 566   mmol/L

    # Figure settings
    fig = plt.figure(figsize=(10, 10))
    
    # Axis settings
    left, bottom, width, height = 0.1, 0.1, 0.8, 0.8
    ax = fig.add_axes([left, bottom, width, height])
    
    # Figure border
    ax.spines['top'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Horizontal lines
    ax.plot([0, 133.4], [0, 0], linestyle='-', linewidth=1.5, color='k')
    ax.plot([0, 133.4], [50, 50], linestyle=':', linewidth=1.5, color='k')
    ax.plot([0, 133.4], [66.7, 66.7], linestyle='-', linewidth=1.25, color='k')
    ax.plot([0, 133.4], [83.4, 83.4], linestyle=':', linewidth=1.5, color='k')
    ax.plot([0, 133.4], [134.4, 134.4], linestyle='-', linewidth=1.5, color='k')
    
    # Vertical lines
    ax.plot([0, 0], [0, 134.4], linestyle='-', linewidth=1.5, color='k')
    ax.plot([50, 50], [0, 134.4], linestyle=':', linewidth=1.5, color='k')
    ax.plot([66.7, 66.7], [0, 134.4], linestyle='-', linewidth=1.25, color='k')
    ax.plot([83.4, 83.4], [0, 134.4], linestyle=':', linewidth=1.5, color='k')
    ax.plot([133.4, 133.4], [0, 134.4], linestyle='-', linewidth=1.5, color='k')
    
    # Limiations
    ax.set_xlim([-20, 145.9])
    ax.set_ylim([-12, 153.4])
    
    # Ticks
    plt.text(0, 135.4, '100', ha='center', va='bottom')
    plt.text(50, 135.4, '50', ha='center', va='bottom')
    plt.text(66.7, 135.4, '33.3', ha='center', va='bottom')
    plt.text(83.4, 135.4, '50', ha='center', va='bottom')
    plt.text(133.4, 135.4, '100', ha='center', va='bottom')
    
    plt.text(-2, 0, '100', ha='right', va='center')
    plt.text(-2, 50, '50', ha='right', va='center')
    plt.text(-2, 66.7, '33.3', ha='right', va='center')
    plt.text(-2, 83.4, '50', ha='right', va='center')
    plt.text(-2, 133.4, '100', ha='right', va='center')
    
    # Lables
    ax.annotate('%' + '$Na^+$' +'+%' + '$K^+$', 
                xy=(0, 145.4), xycoords='data',
                ha="left", va="center",
                xytext=(40, 0), textcoords='offset points',
                size=16,
                arrowprops=dict(arrowstyle="simple", fc="k", ec="k"))
    
    ax.annotate('%' + '$Mg^{2+}$' + '\n'  + '%' + '$Ca^{2+}$', 
                xy=(133.4, 145.4), xycoords='data',
                ha="right", va="center",
                xytext=(-40, 0), textcoords='offset points',
                size=16,
                arrowprops=dict(arrowstyle="simple", fc="k", ec="k"))
    
    ax.annotate('%' + '$Cl^-$', 
                xy=(-12, 0), xycoords='data',
                ha="center", va="bottom",
                xytext=(0, 40), textcoords='offset points',
                size=16, rotation=90,
                arrowprops=dict(arrowstyle="simple", fc="k", ec="k"))
    
    ax.annotate('%' + '$SO4^{2-}$' + '\n' + '%' + '$HCO_3^-$' + '+%' + '$CO_3^{2-}$', 
                xy=(-12, 133.4), xycoords='data',
                ha="center", va="top",
                xytext=(0, -40), textcoords='offset points',
                size=16, rotation=90,
                arrowprops=dict(arrowstyle="simple", fc="k", ec="k"))
    
    # Instrusion arc
    plt.annotate('',
                 xy=(60, 12), xycoords='data',
                 xytext=(180, 200), textcoords='offset points',
                 size=40,
                 # bbox=dict(boxstyle="round", fc="0.8"),
                 arrowprops=dict(arrowstyle="simple",
                                 fc="#F4CED4", ec="none", 
                                 connectionstyle="arc3, rad=-0.3"))
    ax.text(100, 28.4, 'Intrusion', fontsize=28, color="#F4CED4", 
            rotation=45, ha='center', va='center', weight='bold')
             
    # Fresenshing arc
    plt.annotate('',
                 xy=(73.4, 121.4), xycoords='data',
                 xytext=(-180, -200), textcoords='offset points',
                 size=40,
                 # bbox=dict(boxstyle="round", fc="0.8"),
                 arrowprops=dict(arrowstyle="simple",
                                 fc="#D7DFEF", ec="none",
                                 connectionstyle="arc3, rad=-0.3"))
    
    ax.text(33.4, 105, 'Freshening', fontsize=28, color="#D7DFEF", 
            rotation=45, ha='center', va='center', weight='bold')
    
    # Face numbers
    plt.text(25, 108.4, '1', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(25, 75.05, '2', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(25, 58.35, '3', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(25, 25.00, '4', fontsize=26, color="0.6", 
             ha='center', va='center')
    
    plt.text(58.3, 108.4, '5', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(58.3, 75.05, '6', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(58.3, 58.35, '7', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(58.3, 25.00, '8', fontsize=26, color="0.6", 
             ha='center', va='center')
    
    plt.text(75.0, 108.4, '9', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(75.0, 75.05, '10', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(75.0, 58.35, '11', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(75.0, 25.00, '12', fontsize=26, color="0.6", 
             ha='center', va='center')
    
    plt.text(108.4, 108.4, '13', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(108.4, 75.05, '14', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(108.4, 58.35, '15', fontsize=26, color="0.6", 
             ha='center', va='center')
    plt.text(108.4, 25.00, '16', fontsize=26, color="0.6", 
             ha='center', va='center')
    
    
    # Water type notes
    ax.plot([0, 133.4], [-2.5, -2.5], 
            linestyle='-', linewidth=1.0, color='k')
    ax.plot([0, 133.4], [-9.5, -9.5], 
            linestyle='-', linewidth=1.0, color='k')
    ax.plot([50, 50], [-9.5, -2.5], 
            linestyle='-', linewidth=1.0, color='k')
    ax.plot([66.7, 66.7], [-9.5, -2.5], 
            linestyle='-', linewidth=1.0, color='k')
    ax.plot([83.4, 83.4], [-9.5, -2.5], 
            linestyle='-', linewidth=1.0, color='k')
    plt.text(25, -6, '$Na-$', ha='center', va='center', fontsize=13) 
    plt.text(58.3, -6, '$MixNa-$', ha='center', va='center', fontsize=13) 
    plt.text(75, -6, '$MixCa-$', ha='center', va='center', fontsize=13) 
    plt.text(108.4, -6, '$Ca-$', ha='center', va='center', fontsize=13)
    
    ax.plot([135.9, 135.9], [0, 133.4], 
            linestyle='-', linewidth=1.0, color='k')
    ax.plot([142.9, 142.9], [0, 133.4], 
            linestyle='-', linewidth=1.0, color='k')
    ax.plot([135.9, 142.9], [50, 50], 
            linestyle='-', linewidth=1.0, color='k')
    ax.plot([135.9, 142.9], [66.7, 66.7], 
            linestyle='-', linewidth=1.0, color='k')
    ax.plot([135.9, 142.9], [83.4, 83.4], 
            linestyle='-', linewidth=1.0, color='k')
    plt.text(139.4, 108, '$-HCO_3$', 
             ha='center', va='center', fontsize=13, rotation=90) 
    plt.text(139.4, 75.05, '$-MixHCO_3$', 
             ha='center', va='center', fontsize=13, rotation=90) 
    plt.text(139.4, 58.35, '$-MixCl$', 
             ha='center', va='center', fontsize=13, rotation=90) 
    plt.text(139.4, 25, '$-Cl$', 
             ha='center', va='center', fontsize=13, rotation=90) 
    
    
    plt.text(149, 134, 'Hydrochemical Facies', 
             ha='left', va='center', fontsize=16)
    plt.text(149, 127, '1: Na-HCO' + '$_3$' + '/SO' +'$_4$', 
             ha='left', va='center', fontsize=14) 
    plt.text(149, 121, '2: Na-MixHCO' + '$_3$' + '/MixSO' + '$_4$', 
             ha='left', va='center', fontsize=14)
    plt.text(149, 115, '3: Na-MixCl', 
             ha='left', va='center', fontsize=14)
    plt.text(149, 109, '4: Na-Cl', 
             ha='left', va='center', fontsize=14)
    plt.text(149, 103, '5: MixNa-HCO' + '$_3$' + '/SO' +'$_4$', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  97, '6: MixNa-MixHCO' + '$_3$' + '/MixSO' +'$_4$', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  91, '7: MixNa-MixCl', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  85, '8: MixNa-Cl', 
             va='center', fontsize=14)
    plt.text(149,  79, '9: MixCa-HCO' + '$_3$' + '/SO' +'$_4$', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  73, '10: MixCa-MixHCO' + '$_3$' + '/MixSO' +'$_4$', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  67, '11: MixCa-MixCl', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  61, '12: MixCa-Cl', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  55, '13: Ca-HCO' + '$_3$' + '/SO' +'$_4$', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  49, '14: Ca-MixHCO' + '$_3$' + '/MixSO' +'$_4$', ha='left', 
             va='center', fontsize=14)
    plt.text(149,  43, '15: Ca-MixCl', 
             ha='left', va='center', fontsize=14)
    plt.text(149,  37, '16: Ca-Cl', 
             ha='left', va='center', fontsize=14)
    
    # Convert mg/L to meq/L
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
    sumcat = np.sum(meqL[:, 0:4], axis=1)
    suman = np.sum(meqL[:, 4:], axis=1)
    cat = np.zeros((dat.shape[0], 3))
    an = np.zeros((dat.shape[0], 3))
    cat[:, 0] = meqL[:, 0] / sumcat * 100                # Percentage Ca
    cat[:, 1] = meqL[:, 1] / sumcat * 100                # Percentage Mg
    cat[:, 2] = (meqL[:, 2] + meqL[:, 3]) / sumcat * 100 # Percentage Na+K
    an[:, 0] = (meqL[:, 4] + meqL[:, 5]) / suman * 100   # Percentage HCO3 + CO3
    an[:, 2] = meqL[:, 6] / suman * 100                  # Percentage Cl
    an[:, 1] = meqL[:, 7] / suman * 100                  # Percentage SO4
    
    # Convert into cartesian coordinates
    cat_rech = np.where(cat[:, 0] > cat[:, 1], cat[:, 0], cat[:, 1])
    an_rech = np.where(an[:, 0] > an[:, 1], an[:, 0], an[:, 1])
    x =  np.where(cat_rech > cat[:, 2], cat_rech - 100 / 3.0 + 66.7, 100 - cat[:, 2])
    y =  np.where(an_rech > an[:, 2], an_rech - 100 / 3.0 + 66.7, 100 - an[:, 2])
    
    # Plot the scatters
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            plt.scatter(x[i], y[i], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        label=TmpLabel, 
                        edgecolors='black')
        except(ValueError):
                pass
            
    # Calculate the mixing line
    # Coordinates of the left (seawater) of the mixing line
    sumcat_sea = \
        CONC_SWAWATER['Na'] / ions_WEIGHT['Na'] * ions_CHARGE['Na'] + \
            CONC_SWAWATER['K'] / ions_WEIGHT['K'] * ions_CHARGE['K'] + \
                CONC_SWAWATER['Ca'] / ions_WEIGHT['Ca'] * ions_CHARGE['Ca'] + \
                    CONC_SWAWATER['Mg'] / ions_WEIGHT['Mg'] * ions_CHARGE['Mg']
    suman_sea = \
        CONC_SWAWATER['Cl'] / ions_WEIGHT['Cl'] * abs(ions_CHARGE['Cl']) + \
            CONC_SWAWATER['HCO3'] / ions_WEIGHT['HCO3'] * abs(ions_CHARGE['HCO3']) + \
                CONC_SWAWATER['SO4'] / ions_WEIGHT['HCO3'] * abs(ions_CHARGE['HCO3'])          
    x_seawater = \
        100 - (CONC_SWAWATER['Na'] / ions_WEIGHT['Na'] * ions_CHARGE['Na'] + \
            CONC_SWAWATER['K'] / ions_WEIGHT['K'] * ions_CHARGE['K']) / sumcat_sea * 100
    y_seawater = \
        100 - (CONC_SWAWATER['Cl'] / ions_WEIGHT['Cl'] / suman_sea * abs(ions_CHARGE['Cl'])) * 100
    
    # Coordinates of the right (seawater) of the mixing line
    x_rechagrewater = max(x)  # Highest percentage in Ca or Mg
    y_rechagrewater = max(y)  # Highest percentage in HCO3 or SO4
    
    # Plot mixing line
    plt.plot([x_seawater, x_rechagrewater], [y_seawater, y_rechagrewater], 
             linestyle='-', linewidth=1.5, color='k')

    plt.scatter(x_seawater, y_seawater, 
                marker='s', s=75, facecolor='k', edgecolor='k')    
    plt.scatter(x_rechagrewater, y_rechagrewater, 
                marker='s', s=75, facecolor='w', edgecolor='k') 
    
    plt.text(x_seawater + 2.5, y_seawater, 'SW', 
             fontsize=14, ha='left', va='center')
    plt.text(x_rechagrewater + 2.5, y_rechagrewater, 'FW', 
             fontsize=14, ha='left', va='center')

    # Creat the legend
    plt.legend(bbox_to_anchor=(0.15, 0.875), markerscale=1, frameon=False, 
               labelspacing=0.25, handletextpad=0.25)
    
    # Display the info
    print("HFE-D plot created. Saving it now...\n")
    
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
    plot(df, unit='mg/L', figname='HFE-D', figformat='jpg')
    
    
    
    
