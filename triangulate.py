# -*- coding: utf-8 -*-
"""
Created on Mon Jan  6 19:35:55 2025

@author: Eberto Benjumea
"""

import numpy as np
import cv2



P1 = K1 @ np.c_[np.eye(3), np.zeros(3)]
P2 = K2 @ np.c_[R, t]
#print(P1)
#print(P2)
X = cv2.triangulatePoints(P1, P2, pk1, pk2)
X = X[:3]/X[-1] 