# -*- coding: utf-8 -*-
"""
Created on Wed Nov  1 11:43:54 2023

@author: Eberto Benjumea
"""

"""
Created on Sat Oct 28 21:44:08 2023

@author: Eberto Benjumea
"""

import numpy as np
from scipy.interpolate import griddata


# def interpolacion(phase_points_c2, non_nan_coords_c2, sampled_phase_points_c1):
#      corresponding_coords_c2 = griddata(phase_points_c2, non_nan_coords_c2, sampled_phase_points_c1, method='cubic')
#      return corresponding_coords_c2 
res = griddata(points, values, xi, method='cubic')