#====================================================
#          Loading necessary packages
#====================================================

import pyfits
import os
import numpy as np
from astropy.time import Time


#====================================================
#     Loading absolute and temporary variables 
#====================================================

corrpath_data="corrdata.txt"
list_corrdata=list(open(corrpath_data,"r"))
n_corrdata=len(list_corrdata)

errorpath_data="errordata.txt"
list_errordata=list(open(errorpath_data,"r"))


# Hier noch die Positionen und Radien eintragen
target_pos=[2056,1688] #[y,x]
target_rad=[10]
cali1_pos=[2577,1721]
cali1_rad=[8]
cali2_pos=[2036,2288]
cali2_rad=[9]
cali3_pos=[1289,2110]
cali3_rad=[10]

dither=5 #um wie viel pixel das maximum von bild zu bild variieren kann

dcali1_pos=[cali1_pos[0]-target_pos[0], cali1_pos[1]-target_pos[1]]
dcali2_pos=[cali2_pos[0]-target_pos[0], cali2_pos[1]-target_pos[1]]
dcali3_pos=[cali3_pos[0]-target_pos[0], cali3_pos[1]-target_pos[1]]


#=====================================================
#             Definitions for functions
#=====================================================

def dist(x, y, xx, yy):
	x = x * 1.
	y = y * 1.
	xx = xx * 1.
	yy = yy * 1.

	return np.sqrt((yy-y)**2+(xx-x)**2) 


def califlux(tar,cali1,cali2,cali3,etar,ecali1,ecali2,ecali3):
	sigma=((etar/(cali1+cali2+cali3))**2+(tar/(cali1+cali2+cali3)**2)**2*(ecali1**2+ecali2**2+ecali3**2))**0.5
	return sigma


#=====================================================
#                      Main
#=====================================================
path0=list_corrdata[0]
path0=path0.rstrip()
hdulist=pyfits.open(path0) * 1
header0=hdulist[0].header
x_size=header0["NAXIS1"]
y_size=header0["NAXIS2"]

time=[]
date=[]
target_flux=[]
target_error=[]
cali1_flux=[]
cali1_error=[]
cali2_flux=[]
cali2_error=[]
cali3_flux=[]
cali3_error=[]
calibrated_flux=[]
calibrated_error=[]

for i in range(n_corrdata):
	print i

	#target_pos=target_pos #[y,x] Correct Positions for the Image
	#cali1_pos=cali1_pos
	#cali2_pos=cali2_pos
	#cali3_pos=cali3_pos

	path=list_corrdata[i]
	path=path.rstrip()
	hdulist=pyfits.open(path) * 1
	dataset=hdulist[0].data *1
	header=hdulist[0].header

	path2=list_errordata[i]
	path2=path2.rstrip()
	hdulist=pyfits.open(path2) * 1
	dataset_error=hdulist[0].data *1
	header2=hdulist[0].header

	time.append(header["TIME-OBS"])
	date.append(header["DATE-OBS"])


	# Read out Fluxes

	if i==0:
		pos_max=target_pos

	else:
		area=dataset[pos_max[0]-dither:pos_max[0]+dither,pos_max[1]-dither:pos_max[1]+dither] * 1
		new_pos=np.where(area==np.max(area))
		new_pos_a = pos_max[0] - (dither + 1 -new_pos[0]) 
		new_pos_b = pos_max[1] - (dither + 1 -new_pos[1]) 
		pos_max=(new_pos_a,new_pos_b)*1


	print "dataset: ",
	print path
	print "x: ",
	print pos_max[1],
	print "		y: ",
	print pos_max[0]

	target_pos=[pos_max[0] * 1,pos_max[1] * 1] #[y,x] Correct Positions for the Image
	cali1_pos=[target_pos[0] + dcali1_pos[0] ,target_pos[1] + dcali1_pos[1]]
	cali2_pos=[target_pos[0] + dcali2_pos[0] ,target_pos[1] + dcali2_pos[1]]
	cali3_pos=[target_pos[0] + dcali3_pos[0] ,target_pos[1] + dcali3_pos[1]]
	

	#here calc for flux stuff!
	target_index=[]
	for j in range((target_pos[1]-target_rad[0]),(target_pos[1]+target_rad[0]+1)):
		for k in range((target_pos[0]-target_rad[0]),(target_pos[0]+target_rad[0]+1)):
			if dist(target_pos[1],target_pos[0],j,k) <= target_rad[0]:
				target_index.append([k,j])	

	cali1_index=[]
	for j in range((cali1_pos[1]-cali1_rad[0]),(cali1_pos[1]+cali1_rad[0]+1)):
		for k in range((cali1_pos[0]-cali1_rad[0]),(cali1_pos[0]+cali1_rad[0]+1)):
			if dist(cali1_pos[1],cali1_pos[0],j,k) <= cali1_rad[0]:
				cali1_index.append([k,j])	
	cali2_index=[]
	for j in range((cali2_pos[1]-cali2_rad[0]),(cali2_pos[1]+cali2_rad[0]+1)):
		for k in range((cali2_pos[0]-cali2_rad[0]),(cali2_pos[0]+cali2_rad[0]+1)):
			if dist(cali2_pos[1],cali2_pos[0],j,k) <= cali2_rad[0]:
				cali2_index.append([k,j])	

	cali3_index=[]
	for j in range((cali3_pos[1]-cali3_rad[0]),(cali3_pos[1]+cali3_rad[0]+1)):
		for k in range((cali3_pos[0]-cali3_rad[0]),(cali3_pos[0]+cali3_rad[0]+1)):
			if dist(cali3_pos[1],cali3_pos[0],j,k) <= cali3_rad[0]:
				cali3_index.append([k,j])	

	tmp_flux=0.
	tmp_error_flux=0.
	for j in target_index:
		if j[0] >= y_size or j[1] >= x_size:
			tmp_flux=0
			tmp_error_flux=0
			print"ATTENTION: Target or calibrator too close to the limit!!!"
			break
		tmp_flux=tmp_flux + dataset[j[0],j[1]]
		tmp_error_flux=tmp_error_flux + dataset_error[j[0],j[1]]
	target_flux.append(tmp_flux * 1.)
	target_error.append(tmp_error_flux * 1.)

	tmp_flux=0.
	tmp_error_flux=0.
	for j in cali1_index:
		if j[0] >= y_size or j[1] >= x_size:
			tmp_flux=0
			tmp_error_flux=0
			print"ATTENTION: Target or calibrator too close to the limit!!!"
			break
		tmp_flux=tmp_flux + dataset[j[0],j[1]]	
		tmp_error_flux=tmp_error_flux + dataset_error[j[0],j[1]]
	cali1_flux.append(tmp_flux * 1.)
	cali1_error.append(tmp_error_flux * 1.)


	tmp_flux=0.
	tmp_error_flux=0.
	for j in cali2_index:
		if j[0] >= y_size or j[1] >= x_size:
			tmp_flux=0
			tmp_error_flux=0
			print"ATTENTION: Target or calibrator too close to the limit!!!"
			break
		tmp_flux=tmp_flux + dataset[j[0],j[1]]
		tmp_error_flux=tmp_error_flux + dataset_error[j[0],j[1]]	
	cali2_flux.append(tmp_flux * 1.)
	cali2_error.append(tmp_error_flux * 1.)


	tmp_flux=0.
	tmp_error_flux=0.
	for j in cali3_index:
		if j[0] >= y_size or j[1] >= x_size:
			tmp_flux=0
			tmp_error_flux=0
			print"ATTENTION: Target or calibrator too close to the limit!!!"
			break
		tmp_flux=tmp_flux + dataset[j[0],j[1]]	
		tmp_error_flux=tmp_error_flux + dataset_error[j[0],j[1]]
	cali3_flux.append(tmp_flux * 1.)
	cali3_error.append(tmp_error_flux * 1.)
		
	
for i in range(len(target_flux)):
	calibrated_flux.append(target_flux[i]/(cali1_flux[i]+cali2_flux[i]+cali3_flux[i]))	
	calibrated_error.append(califlux(target_flux[i],cali1_flux[i],cali2_flux[i],cali3_flux[i],target_error[i],cali1_error[i],cali2_error[i],cali3_error[i]))

data_output=np.array(zip(time,date,target_flux,target_error,cali1_flux,cali1_error,cali2_flux,cali2_error,cali3_flux,cali3_error,calibrated_flux,calibrated_error), dtype=[("time","S16"), ("date","S16"), ("target_flux",float), ("target_error",float), ("cali1_flux",float), ("cali1_error",float), ("cali2_flux",float), ("cali2_error",float), ("cali3_flux",float), ("cali3_error",float),("calibrated_flux",float), ("calibrated_error",float)])

outname="fluxes2016_04_15.txt"
np.savetxt(outname, data_output, fmt=["%s"]*2 + ["%.5f"]*10, delimiter="|")
print "Your data is saved in '%s'" % (outname)

















