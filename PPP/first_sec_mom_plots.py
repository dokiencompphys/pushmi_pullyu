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
	
	tau1 = 2
	tau2 = tau1
	alpha = [2.1, 1.5, 1.2, 0.6, 0.3, 0.4, 0.2]
	beta = [1.5,2.1, 0.2, 1.7, 1.2, 0.2, 0.5]
	
	N = 8000	# number of realizations for pushmi pullyu paths
	N2 = 10	# number of repeats for msd of msd (error bars)
	
	T = np.arange(0.4, 3, 0.01)
	f = np.zeros(len(T))	
	f2 = np.zeros(len(T))
	
	init=0	# initial state is moving

	c = 2.18

	MSD_s = np.zeros(len(T))
	
	msd_msd = np.zeros(N2)
	
	yerr = np.zeros(len(T))
	
	
	plt.figure()
	plt.title("first moment")
	plt.xlabel('time[min]')
	plt.ylabel('discplacements[um]' )
	
	plt.xscale('log')
	plt.yscale('log')
	for i in np.arange(len(alpha)):
		args = (tau1, tau2, alpha[i], beta[i],c, init)
		f =gav_steh(T, 2, *args)
		plt.plot(T, f, label = r" $\alpha$={:1.2f} and $\beta$={:1.2f}".format(alpha[i], beta[i]), marker='+',markersize = 3.7,  linewidth = 2)
	plt.legend(loc=2, fontsize='large')
	
	plt.figure()
	plt.title("second moment")
	plt.xlabel('time[min]')
	plt.ylabel('discplacements[um]' )
	
	plt.xscale('log')
	plt.yscale('log')
	for i in np.arange(len(alpha)):
		args = (tau1, tau2, alpha[i], beta[i],c, init)
		f =gav_steh(T, 2, *args)
		f2 = gav_steh(T, 3, *args)
		plt.plot(T, f2, label = r" $\alpha$={:1.2f} and $\beta$={:1.2f}".format(alpha[i], beta[i]), marker='+',markersize = 3.7,  linewidth = 2)
	plt.legend(loc=2, fontsize='large')

	plt.show()
	

if __name__=='__main__':
	main()
