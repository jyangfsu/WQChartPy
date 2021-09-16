# -*- coding: utf-8 -*-
"""
Created on Wed Sep 30 20:51:33 2020

@author: Jing
"""
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from .ions import ions_WEIGHT, ions_CHARGE

# Global plot settings
# mpl.rcParams['lines.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 6

def plot_Durvo(df, 
               unit='mg/L', 
               figname='Durvo diagram', 
               figformat='jpg'):
    """Plot the Durvo diagram.
    
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
    # Determine if the required geochemical parameters are defined. 
    if not {'Ca', 'Mg', 'Na', 'K', 'HCO3', 'CO3', 'Cl', 'SO4', 'pH', 'TDS'}.issubset(df.columns):
        raise RuntimeError("""
        Durvo diagram uses geochemical parameters Ca, Mg, Na, K, HCO3, CO3, \
Cl, SO4, pH, and TDS.
        Confirm that these parameters are provided.""")
        
    # Determine if the provided unit is allowed.
    ALLOWED_UNITS = ['mg/L']
    if unit not in ALLOWED_UNITS:
        raise RuntimeError("""
        Currently only mg/L is supported.
        Convert the unit if needed.""")
        
    # Calculate the traingles' location
    h = 0.5 * np.tan(np.pi / 3.0) 
    ltriangle_x = np.array([0, -h, 0, 0])
    ltriangle_y = np.array([0, 0.5, 1, 0])
    ttriangle_x = np.array([0, 0.5, 1, 0])
    ttriangle_y = np.array([1, 1 + h, 1, 1])
    
    # Calculate the rectangles' location 
    crectangle_x = np.array([0, 0, 1, 1, 0]) 
    crectangle_y = np.array([0, 1, 1, 0, 0]) 
    rrectangle_x = np.array([1, 2.618, 2.618, 1, 1]) 
    rrectangle_y = np.array([0, 0, 1, 1, 0]) 
    brectangle_x = np.array([0, 1, 1, 0, 0]) 
    brectangle_y = np.array([0, 0, -0.618, -0.618, 0]) 
    
    # Plot the traingles and rectangles
    fig = plt.figure(figsize=(10,10), dpi=100)
    ax = fig.add_subplot(111, aspect='equal', 
                         frameon=False, xticks=[], yticks=[])
    ax.plot(ltriangle_x, ltriangle_y, '-k', lw=1.0)
    ax.plot(ttriangle_x, ttriangle_y, '-k', lw=1.0)
    # ax.plot(crectangle_x, crectangle_y, '-k', lw=1.0)  # Not needed.
    ax.plot(rrectangle_x, rrectangle_y, '-k', lw=1.0)
    ax.plot(brectangle_x, brectangle_y, '-k', lw=1.0)
    
    # Plot the ticklines with the interval of 0.1 and length of 0.025
    linterval = 0.1
    ticklabels = ['0', '10', '20', '30', '40', '50',
                  '60', '70', '80', '90', '100']
    ticklength = 0.02
    # Left traingle
    for i, x in enumerate(np.linspace(-h, 0, int(1 / linterval + 1))):
        ax.plot([x, x + ticklength * np.sin(np.pi / 3.0)], 
                [-x * np.tan(np.pi / 6.0), -x * np.tan(np.pi / 6.0) + ticklength * np.sin(np.pi / 6.0)], 
                'k', lw=1.0)
        if i in [2, 4, 6, 8]:
            ax.plot([x, 0], [-x * np.tan(np.pi / 6.0), -x / np.sin(np.pi / 3.0)], 'k:')
            ax.text(x - ticklength * np.sin(np.pi / 3.0), -x * np.tan(np.pi / 6.0) - 4 * ticklength * np.sin(np.pi / 6.0), 
                    ticklabels[i], rotation=-30, ha='center', va='center')
    ax.annotate('%' + '$Na^+$' +'+%' + '$K^+$', 
                xy=(-0.2, 0.005), xycoords='data',
                ha="left", va="center",
                xytext=(-100, 32), textcoords='offset points',
                size=12, rotation=-30,
                arrowprops=dict(arrowstyle="simple", fc="k", ec="k"))
    ax.text(-h - 0.05, 0.5, '  100% Mg' + '$^{2+}$', rotation=90, 
            ha='center', va='center', fontsize=12)
    
    for i, x in enumerate(np.linspace(-h, 0, int(1 / linterval + 1))):
        ax.plot([x, x + ticklength * np.sin(np.pi / 3.0)], 
                [0.5 + (h + x) * np.tan(np.pi / 6.0), 0.5 + (h + x) * np.tan(np.pi / 6.0) - ticklength * np.sin(np.pi / 6.0)], 
                'k', lw=1.0)
        if i in [2, 4, 6, 8]:
            ax.plot([x, 0], [0.5 + (h + x) * np.tan(np.pi / 6.0), 1 + x / np.sin(np.pi / 3.0)], 'k:')
            ax.text(x - ticklength * np.sin(np.pi / 3.0), 0.5 + (h + x) * np.tan(np.pi / 6.0) +  4 * ticklength * np.sin(np.pi / 6.0), 
                    ticklabels[i], rotation=30, ha='center', va='center')
    ax.annotate('%' + '$Ca^{2+}$', 
                xy=(-0.3, 1.0), xycoords='data',
                ha="left", va="center",
                xytext=(-68, -28), textcoords='offset points',
                size=12, rotation=30,
                arrowprops=dict(arrowstyle="simple", fc="k", ec="k")) 
    # top traingle
    for i, x in enumerate(np.linspace(0, 0.5, int(1 / linterval + 1))):
        ax.plot([x, x + ticklength * np.sin(np.pi / 6.0)], 
                [1 + x * np.tan(np.pi / 3.0), 1 + x * np.tan(np.pi / 3.0) - ticklength * np.sin(np.pi / 6.0)], 
                'k', lw=1.0)
        if i in [2, 4, 6, 8]:
            ax.plot([x, x / np.sin(np.pi / 6.0)], [1 + x * np.tan(np.pi / 3.0), 1], 'k:')
            ax.text(x - 4 * ticklength * np.sin(np.pi / 6.0), 1 + x * np.tan(np.pi / 3.0) + ticklength * np.sin(np.pi / 6.0),
                    ticklabels[-i - 1], rotation=60, ha='center', va='center')
    ax.annotate('$Cl^-$' + '%', 
                xy=(0.0, 1.28), xycoords='data',
                ha="left", va="center",
                xytext=(15, 42), textcoords='offset points',
                size=12, rotation=60,
                arrowprops=dict(arrowstyle="simple", fc="k", ec="k")) 
    for i, x in enumerate(np.linspace(0.5, 1, int(1 / linterval + 1))):
        ax.plot([x, x - ticklength * np.sin(np.pi / 6.0)], 
                [1 + (1 - x) * np.tan(np.pi / 3.0), 1 + (1 - x) * np.tan(np.pi / 3.0) - ticklength * np.sin(np.pi / 3.0)], 
                'k', lw=1.0)
        if i in [2, 4, 6, 8]:
            ax.plot([x, (x - 0.5) / np.sin(np.pi / 6.0)], [1 + (1 - x) * np.tan(np.pi / 3.0), 1], 'k:')
            ax.text(x +  4 * ticklength * np.sin(np.pi / 6.0), 1 + (1 - x) * np.tan(np.pi / 3.0) + ticklength * np.sin(np.pi / 3.0), 
                    ticklabels[i], rotation=-60, ha='center', va='center')
    ax.annotate('%' + '$HCO_3^-$', 
                xy=(1.00, 1.27), xycoords='data',
                ha="left", va="center",
                xytext=(-49, 49), textcoords='offset points',
                size=12, rotation=-60,
                arrowprops=dict(arrowstyle="simple", fc="k", ec="k")) 
    ax.text(0.5, h + 1.05, '  100% SO' + '$_4^{2-}$', rotation=0, 
            ha='center', va='center', fontsize=12)
    # Center rectangle
    for x in np.linspace(0, 1, int(1 / linterval + 1)):
        ax.plot([x, x], 
                [1, 1 - ticklength], 
                'k', lw=1.0)
    for x in np.linspace(0, 1, int(1 / linterval + 1)):
        ax.plot([0, ticklength], 
                [x, x], 
                'k', lw=1.0)
        
    # Bottom rectangle
    pHlabels = ['6', '6.5', '7', '7.5', '8', '8.5', '9', '9.5']
    for i, x in enumerate(np.linspace(0, -0.618, 8)):
        ax.plot([0, ticklength], 
                [x, x], 
                'k', lw=1.0)
        if i in [2, 4, 6]:
            ax.text(-2 * ticklength, x, pHlabels[i], 
                    ha='center', va='center')
    plt.text(-0.12, -0.618 / 2, 'pH', rotation=90, 
             ha='center', va='center', fontsize=12)
    
    # Right rectangle
    tdslabels = ['0.0', '0.5', '1.0', '1.5', '2.0', '2.5', '3.0', '3.5']
    for i, x in enumerate(np.linspace(1, 2.618, 9)):
        ax.plot([x, x], 
                [0, ticklength], 
                'k', lw=1.0)
        if i in [1, 2, 3, 4, 5, 6, 7]:
            ax.text(x, -2 * ticklength, tdslabels[i], ha='center', va='center')
    ax.text(1 + 1.618 / 2, -0.12, 'TDS (g/L)', ha='center', va='center', fontsize=12)
    
    # Label the watertypes in the central rectangle
    ax.text(0.15, 0.1, '$NaCl$', ha='center', va='center', fontsize=12)
    ax.text(0.8, 0.9, '$CaHCO_3$', ha='center', va='center', fontsize=12)
    
    # Convert chemistry data into plot coordinates
    gmol = np.array([ions_WEIGHT['Ca'], 
                     ions_WEIGHT['Mg'], 
                     ions_WEIGHT['Na'], 
                     ions_WEIGHT['K'], 
                     ions_WEIGHT['HCO3'],
                     ions_WEIGHT['CO3'], 
                     ions_WEIGHT['Cl'], 
                     ions_WEIGHT['SO4'],  
                     1,          # Faked weight for pH
                     1000])      # Faked weight for TDS to convert mg/L to g/L
    eqmol = np.array([ions_CHARGE['Ca'], 
                      ions_CHARGE['Mg'], 
                      ions_CHARGE['Na'], 
                      ions_CHARGE['K'], 
                      ions_CHARGE['HCO3'], 
                      ions_CHARGE['CO3'], 
                      ions_CHARGE['Cl'], 
                      ions_CHARGE['SO4'], 
                      1,         # Faked charge for pH
                      1])        # Faked charge for TDS
    tmpdf = df[['Ca', 'Mg', 'Na', 'K', 'HCO3', 'CO3', 'Cl', 'SO4', 'pH', 'TDS']]
    dat = tmpdf.values

    meqL = (dat / abs(gmol)) * abs(eqmol)
    sumcat = np.sum(meqL[:, 0:4], axis=1)
    suman = np.sum(meqL[:, 4:8], axis=1)
    cat = np.zeros((dat.shape[0], 3))
    an = np.zeros((dat.shape[0], 3))
    cat[:, 0] = meqL[:, 0] / sumcat                     # Percentage Ca
    cat[:, 1] = meqL[:, 1] / sumcat                     # Percentage Mg
    cat[:, 2] = (meqL[:, 2] + meqL[:, 3]) / sumcat      # Percentage Na+K
    an[:, 0] = (meqL[:, 4] + meqL[:, 5]) / suman        # Percentage HCO3 + CO3
    an[:, 2] = meqL[:, 6] / suman                       # Percentage Cl
    an[:, 1] = meqL[:, 7] / suman                       # Percentage SO4
    
    # Convert into cartesian coordinates
    cat_x = -np.sin(np.pi / 3.0) * (1 -  cat[:, 2] - cat[:, 0])
    cat_y = np.sin(np.pi / 6.0) * (1 -  cat[:, 2] - cat[:, 0]) + cat[:, 0]
    an_x = np.sin(np.pi / 6.0) * (1 - an[:, 2]) + np.sin(np.pi / 6.0) * an[:, 0] 
    an_y = 1 + np.sin(np.pi / 3.0) * (1 - an[:, 2] - an[:, 0])
    tds_x = 1 + (meqL[:, 9] - 0) / (4 - 0) * 1.618
    ph_y = -(meqL[:, 8] - 6.0) / (9.5 - 6.0) * 0.618
    
    # Plot the scatters
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            plt.scatter(cat_x[i], cat_y[i], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        #label=TmpLabel, 
                        edgecolors='black')
            plt.scatter(an_x[i], an_y[i], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        #label=TmpLabel, 
                        edgecolors='black')
            plt.scatter(an_x[i], cat_y[i], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        label=TmpLabel, 
                        edgecolors='black')
            plt.scatter(an_x[i], ph_y[i], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        #label=TmpLabel, 
                        edgecolors='black')
            plt.scatter(tds_x[i], cat_y[i], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        #label=TmpLabel, 
                        edgecolors='black')
            
        except(ValueError):
                pass
            
    # Creat the legend
    plt.legend(loc='upper left', markerscale=1, frameon=False, 
               labelspacing=0.25, handletextpad=0.25)
    
    # Save the figure
    plt.savefig(figname + '.' + figformat, format=figformat, 
                bbox_inches='tight', dpi=300)
    
    return

