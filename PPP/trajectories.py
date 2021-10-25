#!usr/bin/python
"""
This program displays
	- trajectories for metastatic cancer cells
	- trajectories for non metastatic cancer cells

"""

import matplotlib.pyplot as plt
import numpy as np
import random
import math
from scipy.special import gamma
from scipy.special import binom
from scipy.optimize import curve_fit
from scipy.stats import moment

from pushmi_pullyu import movement


def main():

	# parameters for the probability distribution
	tau1= 2 # in minutes
	tau2= tau1
	tau3= tau1
	tau4= tau1

	# sub diffusive
	alpha= [0.3, 0.28, 0.35, 0.4, 0.4 ]	#1.7 und 1.5
	beta= [1.8, 1.7, 1.7, 1.8, 2.3]
	gamma= alpha
	delta= beta
	n_sub= len(alpha)
	
	init= 0	# moving as initial state
	
	# super diffusive
	#alpha1= [0.1]
	#beta1= [0.1]
	
	alpha1= [0.7, 0.5, 2.8, 1, 2]
	beta1= [0.3, 0.1, 1,  0.65, 0.4]
	gamma1= alpha1
	delta1= beta1
	n_sup= len(alpha1)

	#velocity	
	c= 2.18 # milli-meter/day
	
	#observation time
	T = 3000

	seeed= 4230
	rng=np.random.default_rng(seeed)
	win_top = 500
	win_bot = -900

	#-------------------------------------------------------------------------------------------------------------
	
	plt.figure(1)
	plt.title("Trajectory of non-metastatic ")
	plt.xlabel('time[min]')
	plt.ylabel('point in space[um]' )
	plt.ylim(bottom=win_bot)
	plt.ylim(top=win_top)
	for i in range(n_sub):
		s1 = random.randint(1,41234)
		s2 = random.randint(1,12311)
		s3 = random.randint(1,9283)
		s4 = random.randint(1,36283)
		plt.plot(movement(alpha[i], beta[i], gamma[i], delta[i], s1, s2, s3, s4, tau1, tau2, tau3, tau4, T, c)[0], movement(alpha[i], beta[i], gamma[i], delta[i], s1, s2, s3, s4, tau1, tau2, tau3, tau4, T, c)[1], marker='.', label =r'$\alpha$={:1.2f}, $\beta$={:1.2f},'.format(alpha[i], beta[i]))
	plt.grid(1)
	plt.legend(loc=3)
	
	plt.show()
	
	rng=np.random.default_rng(seeed)
	
	plt.figure(2)
	plt.title("Trajectory of metastatic ")
	plt.xlabel('time[min]')
	plt.ylabel('point in space[um]' )
	plt.ylim(bottom=win_bot)
	plt.ylim(top=win_top)
	for i in range(n_sup):
		s1 = random.randint(1,41234)
		s2 = random.randint(1,12311)
		s3 = random.randint(1,9283)
		s4 = random.randint(1,36283)
		plt.plot(movement(alpha1[i], beta1[i], gamma1[i], delta1[i], s1, s2, s3, s4, tau1, tau2, tau3, tau4, T, c)[0], movement(alpha1[i], beta1[i], gamma1[i], delta1[i], s1, s2, s3, s4, tau1, tau2, tau3, tau4, T, c)[1], marker='.', label =r'$\alpha$={:1.2f},  $\beta$={:1.2f}'.format(alpha1[i], beta1[i]))
	plt.grid(1)
	plt.legend(loc=3)
	
	plt.show()
	
#-------------------------------------------------------------------------------------------------------------
	
	"""
	plt.figure(3)
	plt.title("Trajectory Random Walker ")
	plt.xlabel('time[min]')
	plt.ylabel('discplacements[um]' )

	i= 0
	s1 = rng.integers(low=10, high=1000, size=1)
	s2 = rng.integers(low=15, high=1000, size=1)
	s3 = rng.integers(low=200, high=1000, size=1)
	s4 = rng.integers(low=439, high=1000, size=1)
	
	time,move1, move2, move3 = movement(alpha[i], beta[i], gamma[i], delta[i], s1, s2, s3, s4, tau1, tau2, tau3, tau4, T, c, init=init)
	
	print(time)
	print(move1)
	print(move2)
	print( move3)
	
	plt.plot(time, move1,\
	  linewidth = 2, color= 'r',marker = '.',  label =r'$\alpha$={:1.2f},  $\beta$={:1.2f}'.format(alpha[i], beta[i]))
	
	plt.plot(time, move2, \
	 linewidth =1.5, linestyle ='--',color= 'b', label ='left leg',marker = '.')
	
	plt.plot(time, move3 , \
	 linewidth = 1.5,linestyle ='--', color='k', label ='right leg', marker = '.')
	
	plt.legend(loc=2)
	
	plt.show()
	"""
if __name__=="__main__":
	main()
	
	
