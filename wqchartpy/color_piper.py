# -*- coding: utf-8 -*-
"""
Created on Wed Oct 22 14:07:24 2021

@author: Luk Peeters
"""
# Load required packages
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pylab import *

from .ions import ions_WEIGHT, ions_CHARGE

# Define the color-coded Piper plotting function
def plot(df, 
         unit='mg/L', 
         figname='Color coded Piper diagram', 
         figformat='jpg'):
    """Plot the color-coded Piper diagram.
    The original color-coded Piper diagram was proposed by Peeters (2014).
    We have updated the background color scheme such that it is more 
    perceptually uniform.
    
    Parameters
    ----------
    df : class:`pandas.DataFrame`
        Geochemical data to draw Gibbs diagram.
    unit : class:`string`
        The unit used in df. Currently only mg/L is supported. 
    color: class `boolean`
        If true, use background coloring of Piper plot.
    alphalevel: class `double`
        Transparency level of points. If 1, points are opaque.
    figname : class:`string`
        A path or file name when saving the figure.
    figformat : class:`string`
        The file format, e.g. 'png', 'pdf', 'svg'
        
    Output
    ----------
    dictionary with:
        cat: [nx3] RGB triple cations
        an:  [nx3] RGB triple anions
        diamond: [nx3] RGB triple central diamond
        
    References
    ----------
    .. [1] Peeters, Luk. 2014.
           A Background Color Scheme for Piper Plots to Spatially Visualize 
           Hydrochemical Patterns
           Groundwater 52(1).
           https://doi.org/10.1111/gwat.12118
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
        
    # Basic shape of piper plot
    # -------------------------------------------------------------------------
    # Define the offset between the diamond and triangle
    alphalevel = 1.0
    offset = 0.10 
    offsety = offset * np.tan(np.pi/3)
    h = 0.5 * np.tan(np.pi/3)
    x1 = 1+2*offset
    
    # Calculate the triangles' location 
    ltriangle_x = np.array([0, 0.5, 1, 0])
    ltriangle_y = np.array([0, h, 0, 0])
    rtriangle_x = ltriangle_x + 2 * offset + 1
    rtriangle_y = ltriangle_y
    
    # Calculate the diamond's location 
    diamond_x = np.array([0.5, 1, 1.5, 1, 0.5]) + offset
    diamond_y = h * (np.array([1, 2, 1, 0, 1])) + (offset * np.tan(np.pi/3))      
    # Array for background
    X  = np.reshape(np.repeat(np.linspace(0, 2 + 2*offset, 1000), 1000 ), (1000, 1000), 'F' )
    Y  = np.reshape(np.repeat(np.linspace(0, 2*h + offsety, 1000), 1000 ), (1000, 1000), 'C' )
    
    # create masks for cation, anion triangle and upper and lower diamond
    ind_cat = np.logical_or(np.logical_and(X<0.5, Y<2*h*X),
                            np.logical_and(X>0.5, Y<(2*h*(1-X))))
    ind_an  = np.logical_or(np.logical_and(X<1.5+(2*offset), Y<2*h*(X-1-2*offset)),
                            np.logical_and(X>1.5+(2*offset), Y<(2*h*(1-(X-1-2*offset)))))
    ind_ld  = np.logical_and(np.logical_or(np.logical_and(X<1.0+offset, Y>-2*h*X + 2*h*(1 + 2*offset)),
                                           np.logical_and(X>1.0+offset, Y>2*h*X - 2*h)),
                             Y < h+offsety)
    ind_ud  = np.logical_and(np.logical_or(np.logical_and( X<1.0+offset, Y <   2*h*X),
                                           np.logical_and( X>1.0+offset, Y <  -2*h*X + 4*h*(1+offset))),
                             Y > h+offsety)
    ind_d   = np.logical_or(ind_ld==1, ind_ud==1)

    # interpolate RGB values
    # rgb_interp = np.load('BivariateColourScheme.npy', allow_pickle=True, fix_imports=True, encoding='latin1').item()
    current_dir = os.path.dirname(__file__)
    rgb_interp = np.load(os.path.join(current_dir, 'BivariateColourScheme.npy'), allow_pickle=True, fix_imports=True, encoding='latin1').item()
    rgba = ['R','G','B','A']
    rgba_dic = {}
    for col in rgba:
        rgba_dic[col] = np.ones((X.shape[0],X.shape[1]))
    # cations
    xc = 2*X[ind_cat]-(Y[ind_cat]/h)-1
    yc = 2*X[ind_cat]+(Y[ind_cat]/h)-1
    for col in rgba[0:3]:
        rgba_dic[col][ind_cat] = np.clip(rgb_interp[col].ev(xc,yc),0,1)
    # anions
    xa = 2*(X[ind_an]-x1) + Y[ind_an]/h - 1
    ya = 2*(X[ind_an]-x1) - Y[ind_an]/h - 1
    for col in rgba[0:3]:
        rgba_dic[col][ind_an] = np.clip(rgb_interp[col].ev(xa,ya),0,1)
    # diamond
    xx = X[ind_d]-(1+offset)
    yy = Y[ind_d]-(offset)
    xd = np.zeros_like(xx)
    yd = np.zeros_like(yy)

    xd[yy>h] = 2*xx[yy>h] -1 + yy[yy>h]/h
    yd[yy>h] = 2*xx[yy>h] - yy[yy>h]/h +1
    xd[yy<=h] = 2*xx[yy<=h] -1 + yy[yy<=h]/h
    yd[yy<=h] = 2*xx[yy<=h] - yy[yy<=h]/h +1
    for col in rgba[0:3]:
        rgba_dic[col][ind_d] = np.clip(rgb_interp[col].ev(xd,yd),0,1)

    RGBA = np.dstack([rgba_dic[a] for a in rgba_dic])

    # Plot the triangles and diamond
    fig = plt.figure(figsize=(10,10), dpi=100)
    ax = fig.add_subplot(111, aspect='equal', frameon=False, 
                         xticks=[], yticks=[])
    ax.imshow(RGBA, 
          origin='lower',
          extent=(0,2+2*offset,0,2*h+offsety))

    ax.plot(ltriangle_x, ltriangle_y, '-k', lw=1.0)
    ax.plot(rtriangle_x, rtriangle_y, '-k', lw=1.0)
    ax.plot(diamond_x, diamond_y, '-k', lw=1.0)
    
    # Plot the dashed lines with the interval of 0.2
    interval = 0.2
    ticklabels = ['0', '20', '40', '60', '80', '100']
    for i, x in enumerate(np.linspace(0, 1, int(1/interval+1))):
        # the left triangle
        ax.plot([x, x - x / 2.0], 
                [0, x / 2.0 * np.tan(np.pi/3)], 
                'k:', lw=1.0)
        ## the bottom ticklabels
        if i in [1, 2, 3, 4]: 
            ax.text(x, 0 - 0.03, ticklabels[-i-1], 
                    ha='center', va='center')
        ax.plot([x, (1 - x) / 2.0 + x], 
                 [0, (1 - x) / 2.0 * np.tan(np.pi/3)], 
                 'k:', lw=1.0)
        ## the right ticklabels
        if i in [1, 2, 3, 4]:
            ax.text((1 - x) / 2.0 + x + 0.026, 
                    (1 - x) / 2.0 * np.tan(np.pi/3) + 0.015, 
                    ticklabels[i], ha='center', va='center', rotation=-60)
        ax.plot([x/2, 1-x/2], 
                [x/2*np.tan(np.pi/3), x/2*np.tan(np.pi/3)], 
                'k:', lw=1.0)
        ## the left ticklabels
        if i in [1, 2, 3, 4]:
            ax.text(x / 2 - 0.026, x / 2 * np.tan(np.pi/3) + 0.015, 
                    ticklabels[i], ha='center', va='center', rotation=60)
        
        # the right triangle
        ax.plot([x + 1 + 2 * offset, x - x / 2.0 + 1 + 2 * offset], 
                [0, x / 2.0 * np.tan(np.pi/3)], 
                'k:', lw=1.0)
        ## the bottom ticklabels
        if i in [1, 2, 3, 4]:
            ax.text(x + 1 + 2 * offset, 0 - 0.03, 
                    ticklabels[i], ha='center', va='center')
        ax.plot([x + 1 + 2 * offset, (1 - x) / 2.0 + x + 1 + 2 * offset],
                 [0, (1 - x) / 2.0 * np.tan(np.pi/3)], 
                 'k:', lw=1.0)
        ## the right ticklabels
        if i in [1, 2, 3, 4]:
            ax.text((1 - x) / 2.0 + x + 1 + 2*offset + 0.026, 
                    (1 - x) /2.0 * np.tan(np.pi/3) + 0.015, 
                    ticklabels[-i-1], ha='center', va='center', rotation=-60)
        ax.plot([x / 2.0 + 1 + 2 * offset, 1 - x / 2 + 1 + 2 * offset], 
                [x / 2.0 * np.tan(np.pi/3), x / 2.0 * np.tan(np.pi/3)], 
                'k:', lw=1.0)
        ## the left ticklabels
        if i in [1, 2, 3, 4]:
            ax.text(x / 2.0 + 1 + 2 * offset - 0.026, 
                    x / 2.0 * np.tan(np.pi/3) + 0.015, 
                    ticklabels[-i-1], ha='center', va='center', rotation=60)
        
        # the diamond
        ax.plot([0.5 + offset + 0.5 / (1/interval) * x/interval, 1 + offset + 0.5 / (1/interval) *x/interval], 
                 [h + offset*np.tan(np.pi/3) + 0.5/(1/interval) * x/interval * np.tan(np.pi/3), offset * np.tan(np.pi/3) + 0.5/(1/interval) * x/interval * np.tan(np.pi/3)], 
                 'k:', lw=1.0)
        ## the upper left and lower right ticklabels
        if i in [1, 2, 3, 4]: 
            ax.text(0.5 + offset +0.5/(1/interval) * x/interval  - 0.026, h + offset * np.tan(np.pi/3) + 0.5/(1/interval) * x/interval * np.tan(np.pi/3) + 0.015, 
                    ticklabels[i], ha='center', va='center', rotation=60)
            ax.text(1 + offset +0.5/(1/interval) * x/interval + 0.026, offset * np.tan(np.pi/3) + 0.5/(1/interval) * x/interval * np.tan(np.pi/3) - 0.015, 
                    ticklabels[-i-1], ha='center', va='center', rotation=60)
        ax.plot([0.5 + offset + 0.5/(1/interval) * x/interval, 1 + offset + 0.5/(1/interval) * x/interval], 
                 [h + offset * np.tan(np.pi/3) - 0.5/(1/interval) * x/interval * np.tan(np.pi/3), 2 * h + offset * np.tan(np.pi/3) - 0.5/(1/interval) * x/interval * np.tan(np.pi/3)], 
                 'k:', lw=1.0)
        ## the lower left and upper right ticklabels
        if i in [1, 2, 3, 4]:  
            ax.text(0.5 + offset + 0.5/(1/interval) * x/interval - 0.026, h + offset * np.tan(np.pi/3) - 0.5/(1/interval) * x/interval * np.tan(np.pi/3) - 0.015, 
                    ticklabels[i], ha='center', va='center', rotation=-60)
            ax.text(1 + offset + 0.5/(1/interval) * x/interval + 0.026, 2 * h + offset * np.tan(np.pi/3) - 0.5/(1/interval) * x/interval * np.tan(np.pi/3) + 0.015, 
                    ticklabels[-i-1], ha='center', va='center', rotation=-60)
    
    # Labels and title
    plt.text(0.5, -offset, '%' + '$Ca^{2+}$', 
             ha='center', va='center', fontsize=12)
    plt.text(1 + 2*offset + 0.5, -offset, '%' + '$Cl^{-}$', 
            ha='center', va='center', fontsize=12)
    plt.text(0.25 - offset*np.cos(np.pi/30), 0.25 * np.tan(np.pi/3) + offset*np.sin(np.pi/30), '%' + '$Mg^{2+}$',  
             ha='center', va='center', rotation=60, fontsize=12)
    plt.text(1.75 + 2*offset + offset*np.cos(np.pi/30), 0.25 * np.tan(np.pi/3) + offset*np.sin(np.pi/30), '%' + '$SO_4^{2-}$',  
              ha='center', va='center', rotation=-60, fontsize=12)
    plt.text(0.75 + offset*np.cos(np.pi/30), 0.25*np.tan(np.pi/3)+offset*np.sin(np.pi/30), '%' + '$Na^+$' + '+%' + '$K^+$',  
              ha='center', va='center', rotation=-60, fontsize=12)
    plt.text(1 + 2*offset + 0.25 - offset*np.cos(np.pi/30), 0.25*np.tan(np.pi/3) + offset*np.sin(np.pi/30), '%' + '$HCO_3^-$' + '+%' + '$CO_3^{2-}$',   
              ha='center', va='center', rotation=60, fontsize=12)
    
    plt.text(0.5 + offset + 0.5*offset + offset*np.cos(np.pi/30), h+offset*np.tan(np.pi/3)+0.25*np.tan(np.pi/3)+offset*np.sin(np.pi/30), '%' + '$SO_4^{2-}$' + '+%' + '$Cl^-$',  
              ha='center', va='center', rotation=60, fontsize=12)
    plt.text(1.5 + offset - 0.25 + offset * np.cos(np.pi/30), h + offset * np.tan(np.pi/3) + 0.25 * np.tan(np.pi/3) + offset*np.sin(np.pi/30), '%' + '$Ca^{2+}$' + '+%' + '$Mg^{2+}$', 
              ha='center', va='center', rotation=-60, fontsize=12)
    
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
    
    # Convert into cartesian coordinates
    cat_x = 0.5 * (2 * cat[:, 2] + cat[:, 1])
    cat_y = h * cat[:, 1]
    an_x = 1 + 2 * offset + 0.5 * (2 * an[:, 2] + an[:, 1])
    an_y = h * an[:, 1]
    d_x = an_y / (4 * h) + 0.5 * an_x - cat_y / (4 * h) + 0.5 * cat_x
    d_y = 0.5 * an_y + h * an_x + 0.5 * cat_y - h * cat_x
    
    # plot data
    plt.plot(cat_x, cat_y, '.k', alpha=alphalevel)
    plt.plot(an_x,  an_y,  '.k', alpha=alphalevel)
    plt.plot(d_x,    d_y,  '.k', alpha=alphalevel)
    
    # calculate RGB values for data
    rgb_dic = {}
    rgb_dic['cat'] = np.zeros((len(cat_x),3))
    xc = 2*cat_x-(cat_y/h)-1
    yc = 2*cat_x+(cat_y/h)-1
    for i,col in enumerate(['R','G','B']):
        rgb_dic['cat'][:,i] = np.clip(rgb_interp[col].ev(xc,yc),0,1)
    rgb_dic['an'] = np.zeros((len(cat_x),3))
    xa = 2*(an_x-x1) + an_y/h - 1
    ya = 2*(an_x-x1) - an_y/h - 1
    for i,col in enumerate(['R','G','B']):
        rgb_dic['an'][:,i] = np.clip(rgb_interp[col].ev(xa,ya),0,1)
    rgb_dic['d'] = np.zeros((len(cat_x),3))
    xx = d_x-(1+offset)
    yy = d_y-(offset)
    xd = np.zeros_like(xx)
    yd = np.zeros_like(yy)
    xd[yy>h] = 2*xx[yy>h] -1 + yy[yy>h]/h
    yd[yy>h] = 2*xx[yy>h] - yy[yy>h]/h +1
    xd[yy<=h] = 2*xx[yy<=h] -1 + yy[yy<=h]/h
    yd[yy<=h] = 2*xx[yy<=h] - yy[yy<=h]/h +1
    for i,col in enumerate(['R','G','B']):
        rgb_dic['d'][:,i] = np.clip(rgb_interp[col].ev(xd,yd),0,1)    
    # Display the info
    print("Color-coded Piper plot created. Saving it now...\n")
    
    # Save the figure
    plt.savefig(figname + '.' + figformat, format=figformat, 
                bbox_inches='tight', dpi=300)
        
    return(rgb_dic)

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
    plot(df, unit='mg/L', figname='color-coded Piper diagram', figformat='jpg')