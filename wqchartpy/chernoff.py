# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 14:33:20 2021

@author: Jing
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *

from .ions import ions_WEIGHT, ions_CHARGE

# Define the Chernoff face plotting function
def plot(df, 
         unit='mg/L', 
         figname='Chernoff face', 
         figformat='jpg'):
    """Plot the Chernoff face.
    
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
    .. [1] Chernoff, H. 1973. 
    The use of faces to represent points in k-dimensional space graphically. 
    Journal of the American Statistical Association 68 no. 342: 361-368. 
    https://doi.org/10.1080/01621459.1973.10482434
    """
    # Basic data check 
    # -------------------------------------------------------------------------
    # Determine if the required geochemical parameters are defined. 
    if not {'Sample', 'Ca', 'Mg', 'Na', 'K', 'HCO3', 'Cl', 'SO4'}.issubset(df.columns):
        raise RuntimeError("""
        Chernoff faces use geochemical parameters Ca, Mg, Na, K, HCO3, Cl, and SO4.
        Also, Sample is requied to save the Chernoff face to local disk for each sample.
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
    
    # Calculate the percentages
    sumcat = np.sum(meqL[:, 0:4], axis=1)
    suman = np.sum(meqL[:, 4:], axis=1)
    cat = np.zeros((dat.shape[0], 3))
    an = np.zeros((dat.shape[0], 3))
    cat[:, 0] = meqL[:, 0] / sumcat                  # Ca
    cat[:, 1] = meqL[:, 1] / sumcat                  # Mg
    cat[:, 2] = (meqL[:, 2] + meqL[:, 3]) / sumcat   # Na+K
    an[:, 0] = meqL[:, 4] / suman                    # HCO3
    an[:, 1] = meqL[:, 5] / suman                    # Cl
    an[:, 2] = meqL[:, 6] / suman                    # SO4
    
    # Plot the Chernoff faces for each sample
    # -------------------------------------------------------------------------
    Labels = []
    for i in range(len(df)):
        if (df.at[i, 'Label'] in Labels or df.at[i, 'Label'] == ''):
            TmpLabel = ''
        else:
            TmpLabel = df.at[i, 'Label']
            Labels.append(TmpLabel)
    
        try:
            x1 = 0.90       # height  of upper face
            x2 = 0.40       # overlap of lower face
            x3 = 0.53       # half of vertical size of face
            
            x4 = cat[i, 1]  # width of upper face, Mg
            x5 = cat[i, 0]  # width of lower face, Ca
            x6 = cat[i, 2]  # length of nose, Na+K
            
            x7 = 0.50       # vertical position of mouth
            x8 = an[i, 2]   # curvature of mouth, SO4
            x9 = an[i, 0]   # width of mouth, HCO3
            
            x10 = 0.73      # vertical position of eyes
            x11 = 0.47      # separation of eyes
            
            x12 = 0.89      # slant of eyes 
            x13 = 0.47      # eccentricity of eyes
            x14 = an[i, 1]  # size of eyes Cl
            x15 = 0.96      # position of pupils
            x16 = 0.98      # vertical position of eyebrows
            x17 = 0.22      # slant of eyebrows
            x18 = 0.27      # size of eyebrows
            
            # 
            fig = plt.figure(figsize=(3,3))
            ax = fig.add_subplot(1,1,1,aspect='equal')
            
            # transform some values so that input between 0,1 yields variety of output
            x3 = 1.9 * (x3 - 0.5)
            x4 = x4 + 0.25
            x5 = x5 + 0.25
            x6 = 0.3 * (x6 + 0.01)
            x8 = 5 * (x8 + 0.001)
            x11 /= 5
            x12 = 2 * (x12 - 0.5)
            x13 += 0.05
            x14 += 0.1
            x15 = 0.5 * (x15 - 0.5)
            x16 = 0.25 * x16
            x17 = 0.5*(x17 - 0.5)
            x18 = 0.5*(x18 + 0.1)

            # Top of face, in box with l=-x4, r=x4, t=x1, b=x3
            e = matplotlib.patches.Ellipse( (0,(x1+x3)/2), 2*x4, (x1-x3), 
                                           fc='white', edgecolor='black', linewidth=2)
            # e.set_clip_box(ax.bbox)
            # e.set_facecolor([0,0,0])
            ax.add_artist(e)
        
            # Bottom of face, in box with l=-x5, r=x5, b=-x1, t=x2+x3
            e = matplotlib.patches.Ellipse( (0,(-x1+x2+x3)/2), 2*x5, (x1+x2+x3), 
                                           fc='white', edgecolor='black', linewidth=2)
            ax.add_artist(e)
        
            # Cover overlaps
            e = matplotlib.patches.Ellipse( (0,(x1+x3)/2), 2*x4, (x1-x3), 
                                           fc='white', edgecolor='black', ec='none')
            ax.add_artist(e)
            e = matplotlib.patches.Ellipse( (0,(-x1+x2+x3)/2), 2*x5, (x1+x2+x3), 
                                           fc='white', edgecolor='black', ec='none')
            ax.add_artist(e)
            
            # Draw nose
            ax.plot([0,0], [-x6/2, x6/2], 'k')
            
            # Draw mouth
            p = matplotlib.patches.Arc( (0,-x7+.5/x8), 1/x8, 1/x8, 
                                       theta1=270-180/pi*arctan(x8*x9), 
                                       theta2=270+180/pi*arctan(x8*x9))
            ax.add_artist(p)
            
            # Draw eyes
            p = matplotlib.patches.Ellipse( (-x11-x14/2,x10), x14, x13*x14, 
                                           angle=-180/pi*x12, 
                                           facecolor='white', edgecolor='black')
            ax.add_artist(p)
            
            p = matplotlib.patches.Ellipse( (x11+x14/2,x10), x14, x13*x14, 
                                           angle=180/pi*x12, 
                                           facecolor='white', edgecolor='black')
            ax.add_artist(p)
        
            # Draw pupils
            p = matplotlib.patches.Ellipse( (-x11-x14/2-x15*x14/2, x10), .05, .05, 
                                           facecolor='black')
            ax.add_artist(p)
            p = matplotlib.patches.Ellipse( (x11+x14/2-x15*x14/2, x10), .05, .05, 
                                           facecolor='black')
            ax.add_artist(p)
            
            # Draw eyebrows
            ax.plot([-x11-x14/2-x14*x18/2,-x11-x14/2+x14*x18/2],
                    [x10+x13*x14*(x16+x17),x10+x13*x14*(x16-x17)],'k')
            ax.plot([x11+x14/2+x14*x18/2,x11+x14/2-x14*x18/2],
                    [x10+x13*x14*(x16+x17),x10+x13*x14*(x16-x17)],'k')
            
            # Show the lables
            ax.text(1.3, 1.2, 'Explanation', ha='left', va='top', fontsize=12)
            ax.text(1.3, 0.9, 'Width of upper face = Mg$^{2+}$', ha='left', va='top', fontsize=12)
            ax.text(1.3, 0.6, 'Width of lower face = Ca$^{2+}$', ha='left', va='top', fontsize=12)
            ax.text(1.3, 0.3, 'Length of nose = Na$^+$+K$^+$', ha='left', va='top', fontsize=12)
            ax.text(1.3, 0.0, 'Curvature of mouth = SO$_4^{2-}$', ha='left', va='top', fontsize=12)
            ax.text(1.3, -0.3, 'Length of mouth = HCO' + '$_3^-$', ha='left', va='top', fontsize=12)
            ax.text(1.3, -0.6, 'Size of eyes = Cl$^-$', ha='left', va='top', fontsize=12)
            
            ax.axis([-1.2, 1.2, -1.2, 1.2])
            ax.set_xticks([])
            ax.set_yticks([])
    
        except(ValueError):
            pass
        
        # Display the info
        print("Chernoff face created for %s. Saving it now...\n" %str(df.at[i, 'Sample']))
    
    
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
    plot(df, unit='mg/L', figname='Chernoff face', figformat='jpg')