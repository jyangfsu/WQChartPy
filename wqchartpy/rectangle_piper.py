# -*- coding: utf-8 -*-
"""
Created on Wed Sep 15 14:54:43 2021

@author: Jing
"""
# Import modules
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *

from .ions import ions_WEIGHT, ions_CHARGE

# Define the plotting function
def plot(df, 
         unit='mg/L', 
         figname='rectaangle Piper diagram', 
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
    .. [1] Ray, R. K. and Mukherjee R. 2008.
           Reproducing the Piper Trilinear Diagram in Rectangular Coordinates.
           Groundwater 46(6), 893–896.
           https://doi.org/10.1111/j.1745-6584.2008.00471.x
    """
    # Basic data check 
    # -------------------------------------------------------------------------
    # Determine if the required geochemical parameters are defined. 
    if not {'Ca', 'Mg', 'Na', 'K', 'HCO3', 'CO3', 'Cl', 'SO4'}.issubset(df.columns):
        raise RuntimeError("""
        Piper diagram uses geochemical parameters Ca, Mg, Na, K, HCO3, CO3, Cl, and SO4.
        Confirm that these parameters are provided.""")
        
    # Determine if the provided unit is allowed
    ALLOWED_UNITS = ['mg/L']
    if unit not in ALLOWED_UNITS:
        raise RuntimeError("""
        Currently only mg/L is supported.
        Convert the unit if needed.""")
        
    # Global plot settings
    # -------------------------------------------------------------------------
    rc('savefig', dpi=300)
    rc('xtick', labelsize=12)
    rc('ytick', labelsize=12)
    rc('font', size=12)
    rc('legend', fontsize=14)
    rc('figure', figsize=(15.0, 5.0)) # define size of Figure window
    markersize=4
    linewidth=2
    xtickpositions = linspace(0, 100, 6) # desired xtickpositions for graphs
    
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

    # Make Figure
    # -------------------------------------------------------------------------
    fig = plt.figure()
            
    # CATIONS
    # -------------------------------------------------------------------------
    ax1 = fig.add_subplot(131)
    ax1b = ax1.twinx()
    
    ax1.fill([100, 0, 100, 100], [0, 100, 100, 0], color = (0.8, 0.8, 0.8))
    ax1.plot([100, 0], [0, 100], 'k')
    ax1.plot([50, 0, 50, 50], [0, 50, 50, 0], 'k--')
    ax1.text(25, 15, 'Na type', fontsize=14)
    ax1.text(75, 15, 'Ca type', fontsize=14)
    ax1.text(25, 65, 'Mg type', fontsize=14)
    
    # Plot the scatters
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            plt.scatter(100 * cat[i, 0], 100 * cat[i, 1], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        label=TmpLabel, 
                        edgecolors='black')
        except(ValueError):
                pass
            
    # Creat the legend
    lgnd = legend(loc='upper left', markerscale=1, frameon=False, fontsize=12,
                  handletextpad=-0.5, facecolor=(0.8, 0.8, 0.8))
    #lgnd.legendHandles[0]._sizes = [30]
    #lgnd.legendHandles[1]._sizes = [30]
    
    # Set first axis limits and labels
    ax1.set_xlim(0,100)
    ax1.set_ylim(0,100)
    ax1b.set_ylim(0,100)
    ax1.set_xlabel('$\longleftarrow$' + 'Ca$^{2+}$ (% meq)', 
                   fontsize=12, weight='normal')
    ax1b.set_ylabel('Mg$^{2+}$ (% meq)' + '$\longrightarrow$', 
                    fontsize=12, weight='normal')
    #ax1.set_xlabel('<- Ca (% meq)')
    #ax1b.set_ylabel('Mg (% meq) ->')
    setp(ax1, yticklabels=[])
    
    # Reverse x axis:
    ax1.set_xlim(ax1.get_xlim()[::-1]) 
    
    # ANIONS
    # -------------------------------------------------------------------------
    ax3 = fig.add_subplot(133)
    ax3.fill([100, 100, 0, 100], [0, 100, 100, 0], color = (0.8, 0.8, 0.8))
    ax3.plot([0, 100], [100, 0], 'k')
    ax3.plot([50, 50, 0, 50], [0, 50, 50, 0],'k--')
    ax3.text(55, 15, 'Cl type', fontsize=14)
    ax3. text(5, 15, 'HCO$_3$ type', fontsize=14)
    ax3.text(5, 65, 'SO$_4$ type', fontsize=14)
    
    # Plot the scatters
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            plt.scatter(100 * an[i, 2], 100 * an[i, 1], 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        #label=TmpLabel, 
                        edgecolors='black')
        except(ValueError):
                pass
            
    # Set first axis limits and labels        
    ax3.set_xlim(0, 100)
    ax3.set_ylim(0, 100)
    ax3.set_xlabel('Cl$^-$ (% meq)' + '$\longrightarrow$', 
                   size=12, weight='normal')
    ax3.set_ylabel('SO$_4^{2-}$ (% meq)' + '$\longrightarrow$', 
                   size=12, weight='normal')
    #ax3.set_xlabel('Cl (% meq) ->')
    #ax3.set_ylabel('SO4 (% meq) ->')
    
    # CATIONS AND ANIONS COMBINED IN DIAMOND SHAPE PLOT
    # -------------------------------------------------------------------------
    ax2 = fig.add_subplot(132)
    ax2b = ax2.twinx()
    
    # Create plot
    ax2.plot([0, 100], [10, 10], 'k--')
    ax2.plot([0, 100], [50, 50], 'k--')
    ax2.plot([0, 100], [90, 90], 'k--')
    ax2.plot([10, 10], [0, 100], 'k--')
    ax2.plot([50, 50], [0, 100], 'k--')
    ax2.plot([90, 90], [0, 100], 'k--')
    
    # Plot the scatters
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            ax2.scatter(100 * cat[i, 2], 100 * (an[i, 1] + an[i, 2]), 
                        marker=df.at[i, 'Marker'],
                        s=df.at[i, 'Size'], 
                        color=df.at[i, 'Color'], 
                        alpha=df.at[i, 'Alpha'],
                        #label=TmpLabel, 
                        edgecolors='black')
        except(ValueError):
                pass
            
    # Set second axis limits and labels
    ax2.set_xlim(0,100)
    ax2.set_ylim(0,100)
    ax2.set_xlabel('Na$^+$+K$^+$ (% meq)' + '$\longrightarrow$', 
                   size=12, weight='normal')
    ax2.set_ylabel('SO$_4^{2-}$+Cl$^-$ (% meq)' + '$\longrightarrow$', 
                   size=12, weight='normal')
    ax2.set_title('$\longleftarrow$' + 'Ca$^{2+}$+Mg$^{2+}$ (% meq)', 
                  size=12, weight='normal')
    ax2b.set_ylabel('$\longleftarrow$' + 'CO$_3^{2-}$+HCO$_3^-$ (% meq)', 
                    size=12, weight='normal')
    
    #ax2.set_xlabel('Na+K (% meq) ->')
    #ax2.set_ylabel('SO4+Cl (% meq) ->')
    #ax2.set_title('<- Ca+Mg (% meq)', fontsize = 12)
    #ax2b.set_ylabel('<- CO3+HCO3 (% meq)')
    #ax2b.set_ylim(0,100)
    
    # Reverse 2nd y axis:
    ax2b.set_ylim(ax2b.get_ylim()[::-1])
    
    # adjust position of subplots
    plt.subplots_adjust(left=0.05, bottom=0.2, right=0.95, top=0.90, 
                    wspace=0.4, hspace=0.0)
    
    # Display the info
    print("Rectangle Piper plot created. Saving it now...\n")
    
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
    plot(df, unit='mg/L', figname='rectangle Piper diagram', figformat='jpg')
    
    
    
    
    
    
    
    
    
    
    
    
   