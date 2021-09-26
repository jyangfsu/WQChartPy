# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 13:41:09 2021

@author: Jing
"""
# Weight values are taken from hanford.dat provided by PFLOTRAN 
#   https://pflotran.org/.

ions_WEIGHT = {'Ca'  : 40.0780,
               'Mg'  : 24.3050,
               'K'   : 39.0983,
               'Na'  : 22.9898,
               'Cl'  : 35.4527,
               'SO4' : 96.0636,
               'CO3' : 60.0092,
               'HCO3': 61.0171}

ions_CHARGE = {'Ca'  : +2,
               'Mg'  : +2,
               'K'   : +1, 
               'Na'  : +1,
               'Cl'  : -1,
               'SO4' : -2,
               'CO3' : -2,
               'HCO3': -1,}
    