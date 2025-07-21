# -*- coding: utf-8 -*-
# """
# Created on Mon Nov 27 10:13:23 2023

# @author: Eberto Benjumea
# """

import numpy as np
import cv2
#import matplotlib.pyplot as plt


# Load the .mat file (Replace 'data_correspondence.mat' with the path to your .mat file)
K1
K1_dist
M
worldPoints
correct_poses
res_K1
res_K2

#plt.figure();
#plt.scatter(worldPoints[:,0],worldPoints[:,1],marker="o", color='blue')
#plt.show()

#%%
print(M.shape)
print(correct_poses.shape)
correct_poses=correct_poses-1

correct_poses=correct_poses.ravel()

correct_poses = correct_poses.tolist()
correct_poses = [int(valor) for valor in correct_poses]
#print(correct_poses)
M_K1=M[:,:,:,0]
M_K2=M[:,:,:,1]
#print(M_K1.shape)
#print(M_K2.shape)

M_K1=M_K1[:,:,correct_poses]
M_K2=M_K2[:,:,correct_poses]
print(correct_poses)
print(M_K1.shape)
print(M_K2.shape)

worldPoints=np.float32(worldPoints)
cero=np.expand_dims(np.zeros(worldPoints.shape[0]), axis=0)
worldPoints=np.concatenate((worldPoints,cero.T), axis=1)
M_K1=np.float32(M_K1)
M_K2=np.float32(M_K2)


# view_points=0
# if view_points==1:
#     for i in range(0,M_K1.shape[2]):
#         fig, axes = plt.subplots(1, 2, figsize=(12, 6)) 
#         axes[0].imshow(np.zeros((res_K1[0],res_K1[1])))
#         axes[0].plot(M_K1[:,0,i],M_K1[:,1,i], 'ro', markersize=2)
#         axes[0].set_title('K1 points')
#         axes[0].axis('off')
#         axes[1].imshow(np.zeros((res_K2[0],res_K2[1])))
#         axes[1].plot(M_K2[:,0,i],M_K2[:,1,i], 'ro', markersize=2)
#         axes[1].set_title('K2 points')
#         axes[1].axis('off')
#         #plt.show()
#         plt.pause(1)
#         plt.close('all')
  
M_K1l=list()
M_K2l=list() 
worldPointsl=[] 

for i in range(0,M_K1.shape[2]):
    M_K1l.append(M_K1[:,:,i].tolist())
    M_K2l.append(M_K2[:,:,i].tolist())
    worldPointsl.append(worldPoints.tolist())

# for i in range(0,M_K1.shape[2]):
#     M_K1l.append(M_K1[i])
#     M_projl.append(M_proj[i].T)
#     worldPointsl.append(worldPoints)

objpoints=np.array(worldPointsl, dtype = np.float32)
K1points=np.array(M_K1l, dtype = np.float32)
K2points=np.array(M_K2l, dtype = np.float32)


CALIBFLAG = 0  # cv2.CALIB_FIX_K3
#print(' ', len(objp_list), 'images are used')
# rms1, K1, K1_dist, rvecs, tvecs = cv2.calibrateCamera(
#     objpoints , K1points, (int(res_K1[0]),int(res_K1[1])), None, None, flags = cv2.CALIB_FIX_K3)

rms2, K2, K2_dist, rvecs, tvecs = cv2.calibrateCamera(
    objpoints , K2points, (int(res_K2[0]),int(res_K2[1])), None, None, flags = cv2.CALIB_FIX_K3)
print('Bien')
#%%
k1copy=K1.copy()
k2copy=K2.copy()
ret, K1, D1, K2, D1, R, T, E, F=cv2.stereoCalibrate(objpoints, K1points, K2points, k1copy, K1_dist, k2copy, K2_dist, None,flags=cv2.CALIB_FIX_INTRINSIC+cv2.CALIB_FIX_K3)
print('Bien')
data = {"K1": K1, "K1_dist":K1_dist, "K2":K2, "K2_dist":K2_dist, "R":R, "t":T,"E":E, "F":F,"rms2":rms2, "rett":ret}
