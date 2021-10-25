#!usr/bin/python
"""
This program computes the
	- the inverse laplace transformation of certain functions with gaver stehfest

"""
from scipy.special import gammaincc
from scipy.special import binom
from scipy.special import gamma
from scipy.optimize import curve_fit
import random
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from gaver_steh import sec_mom_asy, gav_steh

from pushmi_pullyu import create_hist
	

def main():
	
	tau1 = 2
	tau2 = tau1
	alpha = 0.4
	beta = 0.2
	
	N = 5000	# number of realizations for pushmi pullyu paths
	N2 = 10

	T = np.arange(10, 3000, 500	)
	f = np.zeros(len(T))	
	f2 = np.zeros(len(T))

	c = 2.18

	MSD_s = np.zeros(len(T))

	msd_msd = np.zeros(N2)
	
	yerr = np.zeros(len(T))

	print("T		mean_s		MSD_s \n")
	print("-------------------------------------\n")
	for k in np.arange(len(T)):
		# imported functions to calculate MSD and mean for one leg
		MSD_s[k] = create_hist(alpha, beta, alpha, beta, tau1, tau2, tau1, tau2, T[k], c, N)[4]
		print("{:3.2f}	{:3.2f}	 \n".format(T[k], MSD_s[k] ))

		for l in np.arange(N2):
			 msd_msd[l]=create_hist(alpha, beta, alpha, beta, tau1, tau2, tau1, tau2, T[k], c, N)[4]
		yerr[k] = np.sqrt(np.var(msd_msd))

	plt.figure()
	plt.title(r"$ \langle x^2 \rangle (t) $ with $\alpha$={:1.2f} and $\beta$={:1.2f}".format(alpha, beta))
	plt.xlabel('time[min]')
	plt.ylabel('discplacements[um]' )
	
	plt.xscale('log')
	plt.yscale('log')

	plt.errorbar(T, MSD_s, label='simulation model', marker = '*',yerr =yerr)
	plt.grid('True')
	plt.legend(loc=2, fontsize='large')

	plt.show()

if __name__=='__main__':
	main()
