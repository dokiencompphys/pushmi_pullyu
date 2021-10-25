#!usr/bin/python
"""
This program tests for low values alpha wether the gaver stehfest of the pdf in laplace domain matches the pdf in time domain

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
	
	
def psi(exp, tau, t):
	return exp * tau**exp / (tau+t)**(1+exp)

def main():
	
	tau1 = 2
	tau2 = tau1
	alpha = [2.1, 1.5, 1.2, 0.6, 0.1, 0.4, 0.2]
	beta = [1.5,2.1, 0.2, 1.7, 1.2, 0.2, 0.5]

	
	T = np.arange(10, 3200, 10)
	f = np.zeros(len(T))	
	f2 = np.zeros(len(T))

	c = 2.18

	MSD_s = np.zeros(len(T))

	for i in np.arange(len(alpha)):
		args = (tau1, tau2, alpha[i], beta[i],c)
		f2 = gav_steh(T, 1, *args)			#pdf in laplace

		plt.figure()
		plt.title(r"$ \langle x^2 \rangle (t) $ with $\alpha$={:1.2f} and $\beta$={:1.2f}".format(alpha[i], beta[i]))
		plt.xlabel('time[min]')
		plt.ylabel('discplacements[um]' )
		
		plt.errorbar(T, f2, label = 'gaver - stehfest on pdf in laplace', marker='+',markersize = 6.7,  linewidth = 2)
		plt.errorbar(T, psi(alpha[i],tau1, T), label='pdf in time domain', marker ='o', markersize = 5, linewidth = 2)
		plt.legend(loc=2, fontsize='large')

	plt.show()

if __name__=='__main__':
	main()
