#!usr/bin/python
"""
This program computes the
	- the inverse laplace transformation of certain functions with gaver stehfest

"""

import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from gaver_steh import sec_mom_asy, gav_steh

from pushmi_pullyu import create_hist

import time
	

def main():
	
	start = time.time()

	tau1 = 2
	tau2 = tau1
	alpha = [2.1, 1.5, 1.2, 0.6, 0.3, 0.4, 0.2]
	beta = [1.5,2.1, 0.2, 1.7, 1.2, 0.2, 0.5]
	
	N = 8000	# number of realizations for pushmi pullyu paths
	N2 = 10	# number of repeats for msd of msd (error bars)
	
	T = np.arange(10, 11000, 1000)
	f = np.zeros(len(T))	
	f2 = np.zeros(len(T))
	
	init=0	# initial state is moving

	c = 2.18

	MSD_s = np.zeros(len(T))
	
	msd_msd = np.zeros(N2)
	
	yerr = np.zeros(len(T))

	#for i in np.arange(len(alpha)):
	for i in np.arange(5,6):
		args = (tau1, tau2, alpha[i], beta[i],c, init)
		f =gav_steh(T, 2, *args)
		f2 = gav_steh(T, 3, *args)

		asymp2mom = sec_mom_asy(2, T, *args)

		print("T		MSD_gav		MSD_num	yerr \n")
		print("-------------------------------------\n")
		for k in np.arange(len(T)):
			
			
			for l in np.arange(N2):
				msd_msd[l]=create_hist(alpha[i], beta[i], alpha[i], beta[i], tau1, tau2, tau1, tau2, T[k], c, N, init =init)[1]
			
			yerr[k] = np.sqrt(np.var(msd_msd))
			
			# imported functions to calculate MSD and mean for one leg
			MSD_s[k] = np.mean(msd_msd)
			print("{:3.2f}	{:5.2f}	{:5.2f}	{:3.2f}	 \n".format(T[k], f2[k], MSD_s[k], yerr[k] ))
	
		plt.figure()
		plt.title(r"$ \langle x^2 \rangle (t) $ with $\alpha$={:1.2f} and $\beta$={:1.2f}".format(alpha[i], beta[i]))
		plt.xlabel('time[min]')
		plt.ylabel('discplacements[um]' )
		
		plt.xscale('log')
		plt.yscale('log')

		plt.errorbar(T, f2, label = 'gaver - stehfest', marker='+',markersize = 6.7,  linewidth = 2)
		plt.errorbar(T, asymp2mom, label='asymptotics', marker ='o', markersize = 5, linewidth = 2)
		plt.errorbar(T, MSD_s, label='simulation model', marker = '*', markersize = 6.7,yerr = yerr, fmt = 'o', capsize=1)
		plt.legend(loc=2, fontsize='large')

	plt.show()
	
	end = time.time()
	print("The program took {:3.2f}".format(end-start))

if __name__=='__main__':
	main()
