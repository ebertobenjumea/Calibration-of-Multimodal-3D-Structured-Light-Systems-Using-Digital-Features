# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 10:13:23 2023

@author: Eberto Benjumea
"""

import numpy as np
import cv2
import matplotlib.pyplot as plt
from scipy.io import loadmat
from scipy.io import savemat


# Load the .mat file (Replace 'data_correspondence.mat' with the path to your .mat file)
mat_contents = loadmat('data_calib_two_stages_calib2.mat')


# Extract relevant variables
M= mat_contents['M']
worldPoints = mat_contents['worldPoints']

plt.figure();
plt.scatter(worldPoints[:,0],worldPoints[:,1],marker="o", color='blue')
plt.show()

print(M.shape)
M_vis=M[:,:,:,0]
M_ir=M[:,:,:,1]

width=2048
height=1600

worldPoints=np.float32(worldPoints)
cero=np.expand_dims(np.zeros(worldPoints.shape[0]), axis=0)
worldPoints=np.concatenate((worldPoints,cero.T), axis=1)
M_vis=np.float32(M_vis)
M_ir=np.float32(M_ir)


for i in range(0,M_vis.shape[2]):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6)) 
    axes[0].imshow(np.zeros((height,width)))
    axes[0].plot(M_vis[:,0,i],M_vis[:,1,i], 'ro', markersize=2)
    axes[0].set_title('VIS points')
    axes[0].axis('off')
    axes[1].imshow(np.zeros((192,256)))
    axes[1].plot(M_ir[:,0,i],M_ir[:,1,i], 'ro', markersize=2)
    axes[1].set_title('IR points')
    axes[1].axis('off')
    #plt.show()
    plt.pause(1)
    plt.close('all')
    
M_visl=list()
M_irl=list() 
worldPointsl=[]    
for i in range(0,M_vis.shape[2]):
    if (sum(sum(M_vis[:,:,i]))==0.0 or  sum(sum(M_ir[:,:,i]))==0.0)==False:
        M_visl.append(M_vis[:,:,i].tolist())
        M_irl.append(M_ir[:,:,i].tolist())
        worldPointsl.append(worldPoints.tolist())

# for i in range(0,M_vis.shape[2]):
#     M_visl.append(M_vis[i])
#     M_irl.append(M_ir[i].T)
#     worldPointsl.append(worldPoints)

objpoints=np.array(worldPointsl, dtype = np.float32)
imgpoints=np.array(M_visl, dtype = np.float32)
irpoints=np.array(M_irl, dtype = np.float32)


# ---------------------Cleaning data--------------------------- 
# cleaned_imgpoints=list()
# cleaned_irpoints=list()
# cnt_useful_poses=0
# for i in range(0,imgpoints.shape[0]):
#     if (sum(sum(imgpoints[i]))==0.0 or  sum(sum(irpoints[i]))==0.0)==False:
#         #print(i)
#         #print(cnt_useful_poses)
#         cleaned_imgpoints.append(imgpoints[i])
#         cleaned_irpoints.append(irpoints[i])
#         cnt_useful_poses += 1
#     else:
#         objpoints.pop(i)
        
# imgpoints= cleaned_imgpoints   
# irpoints=cleaned_irpoints    

#--------------------------------------------------------------    
# Load the .mat file with intrinsic marices(Replace 'data_correspondence.mat' with the path to your .mat file)
mat_contents2 = loadmat('intrinsic_parameters.mat')


# Extract relevant variables
k1= mat_contents2['K1']
kt= mat_contents2['Kt']
dist1= mat_contents2['dist1']
dist12= mat_contents2['dist12']
distir_1= mat_contents2['distir_1']
distir_2= mat_contents2['distir_2']

cam_dist=np.concatenate((dist1,dist12), axis=1)
t_dist=np.concatenate((distir_1,distir_2), axis=1)

#----------------------------------------------------------

CALIBFLAG = 0  # cv2.CALIB_FIX_K3
vis_shape=[height, width]
#print(' ', len(objp_list), 'images are used')
rms, k1, cam_dist, rvecs, tvecs = cv2.calibrateCamera(
   objpoints , imgpoints, (height, width), None, None)

rms, kt, t_dist, rvecs, tvecs = cv2.calibrateCamera(
   objpoints , irpoints, (192, 256), None, None)


ret, K1, D1, Kt, Dt, R, T, E, F=cv2.stereoCalibrate(objpoints, imgpoints, irpoints, k1, cam_dist, kt, t_dist, None,flags=cv2.CALIB_FIX_INTRINSIC+cv2.CALIB_FIX_K2)


# Save data in .mat file
data = {"k1": k1, "cam_dist": cam_dist, "kt": kt, "t_dist": t_dist,"R": R, "T": T}


savemat("data_from_python_calib2.mat", data)

