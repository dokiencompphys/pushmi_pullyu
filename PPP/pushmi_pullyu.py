#!usr/bin/python
"""
This program provides the functions 
	- to draw a wating time from a PDF
	- combine waiting time strings from two until an observation time
	- create trajectory of the pushmi-pullyu-animal

"""

import numpy as np
import random


def draw_wait_time_pareto(x, exp, tau):
	"""
	returns waiting time for a pareto 2 distribution.
	"""
	return tau/(1-x)**(1/exp) - tau

def find_n(T, exp, tau):
	"""
	return index n. first generate the random number sequence and find out
	how many elements n it needs for the waiting time sum to toal time T.
	"""
	T_i= 0
	n = 0
	while T_i < T:
		x=np.random.random_sample()
		t = draw_wait_time_pareto(x,exp,tau)
		T_i= T_i + t
		n = n + 1 
	return n

def get_memory_string(exp1, exp2,s1, s2, tau1, tau2, T, init=0):
	"""
	return memory string of waiting time history,
	sum totals over T
	returns time of total movement up to observation T
	"""
	#draw waiting times until exceeds T
	#this step is to ensure we have enough waiting times drawn
	np.random.seed(s1)
	n1 = find_n(T, exp1, tau1)
	
	np.random.seed(s2)	
	n2 = find_n(T, exp2, tau2)
	
	if init ==0:
		#reset seeds to ensure the same waiting time sequences like in find_n()
		# and store waiting times
		np.random.seed(s1)
		x1 = np.random.random_sample(n1)
		Wt1 = draw_wait_time_pareto(x1, exp1, tau1)
		
		np.random.seed(s2)
		x2 = np.random.random_sample(n2)
		Wt2 = draw_wait_time_pareto(x2, exp2, tau2)
	elif init ==1:
	# swap odd with even ordered élements in array
		np.random.seed(s1)
		x1 = np.random.random_sample(n1)
		Wt2 = draw_wait_time_pareto(x1, exp1, tau1)
		
		np.random.seed(s2)
		x2 = np.random.random_sample(n2)
		Wt1 = draw_wait_time_pareto(x2, exp2, tau2)
	
	#combine both waiting times
	Wt = np.zeros (n1+n2)
	
	T_i=0
	i = 0
	t_1 = 0
	t_2 = 0
	while T_i < T:				# if equal, then it will stop too
		Wt[2*i ] = Wt1[i]
		Wt[2*i + 1 ] = Wt2[i]

		T_i = T_i + Wt1[i] + Wt2[i]	# total current time
		t_1 = t_1 + Wt1[i]
		t_2 = t_2 + Wt2[i]
		 
		i = i + 1

	i = i -1	# reduce to get the correct last index
		
	# cut elements if exceed total time T
	# three cases can happen, end of time falls in state 1 or 2
	
	# T falls into state 2, because the last element in waiting time array is needed
	# to exceed observation time T
	if T_i - Wt[2*i+1] < T:			
		t_1 = sum(Wt1[:i+1])
		t_2 = T - t_1	
		Wt[2*i+1] = t_2 - sum(Wt2[:i])
		t_move = t_2
		t_rest = t_1																				
		#print("state 2")
		
		if init == 0:
			t_move = t_2
			t_rest = t_1
		elif init == 1:
			t_move = t_1
			t_rest = t_2
		
		return Wt[:2*i+2], t_move, t_rest + t_move, T, sum(Wt[:2*i+2]), len(Wt[:2*i])
		# to access the desired range in an array in python we need to add one to the index 
	
	# T falls into state 1	
	elif T_i - Wt[2*i+1] > T:
		t_2 = sum(Wt2[:i])
		t_1 = T- t_2 
		Wt[2*i] = t_1 - sum(Wt1[:i])
		#print("state 1")
		
		if init == 0:
			t_move = t_2
			t_rest = t_1
		elif init == 1:
			t_move = t_1
			t_rest = t_2
		
		return Wt[:(2*i)+1], t_move, t_rest + t_move, T, sum(Wt[:2*i+1]), len(Wt[:2*i+1])
	
	# falls right between state 1 and 2
	elif T_i - Wt[2*i+1] == T:
		t_2 = t_2 - Wt[2*i+1]
		if init == 0:
			t_move = t_2
			t_rest = t_1
		elif init == 1:
			t_move = t_1
			t_rest = t_2
		
		return Wt[:(2*i)+1], t_move, t_rest + t_move, T, sum(Wt[:2*i+1]), len(Wt[:2*i+1])
	

def movement(exp1, exp2, exp3, exp4, s1, s2, s3, s4, tau1, tau2, tau3, tau4, T, c, init=0):
	"""
	
	"""
	Wt_left = get_memory_string(exp1, exp2,s1, s2, tau1, tau2, T, init = init)[0]
	Wt_right = get_memory_string(exp3, exp4,s3, s4, tau3, tau4, T, init= init)[0] 
	
	print(Wt_left)
	print(Wt_right)
	
	N = len(Wt_left) + 1 # the waiting times are M+N edges, space points therefor need to be M+N+1 
	M = len(Wt_right) + 1
	Tx = np.zeros(M + N)
	X = np.zeros(M + N) # total point of space

	#for single legs
	Xl = np.zeros(M + N)
	Xr = np.zeros(M + N)
	
	#index of waiting time array
	l =0
	r= 0

	timer_l = Wt_left[0]
	timer_r = Wt_right[0]

	if init==2:
		init = random.randint(0,1)	# randomized move or resting state

	if timer_l > timer_r:
	# that means the right leg will produce the first time point
		l0=init
		r0=init+1		# in the latter while loop 	
					# r will be increased, we correct it by 
					# introducing this phase shift
	elif timer_l < timer_r:
		l0=init+1
		r0=init
	else:
		l0=init+1
		r0=init+1
		
	n=0
	#introduce timers for left and right leg
	while T-Tx[n] > 10**(-9):
		if timer_l - timer_r > 0:
			#if the remaining waiting time on 
			#the left leg is larger, add the waiting time
			#of the right string to the total observation time

			Tx[n+1] = Tx[n] +  timer_r 	#next increment in running time
			timer_l = timer_l - timer_r
			Dt = timer_r			# determines the stretch of PPM

			r = r + 1 #next data point for right memory string
			if r < M-1:
				timer_r = Wt_right[r] # take next waiting time 
			

		elif timer_l - timer_r < 0:
			
			Tx[n+1] = Tx[n] + timer_l
			timer_r = timer_r - timer_l
			Dt = timer_l
			l = l + 1	
			if l < N-1:
				timer_l = Wt_left[l]
				
		else:
		# timer_l = timer_r
			Tx[n+1] = Tx[n] + timer_l
			Dt = timer_l
			l = l + 1
			r = r+1
			if l < N-1:
				timer_l = Wt_left[l]
			if r < M-1:
				timer_r = Wt_right[r]

		Xl[n+1] = Xl[n] - (l+l0)%2 *c * Dt
		Xr[n+1] = Xr[n] + (r+r0)%2 *c * Dt
		X[n+1] = Xr[n+1] + Xl[n+1]

		n = n+1
		#print("n={:1.1f} Tx={:f} Xl={:f} Xr = {:f} X ={:f} l={:f} r={:f} \n".format(n, Tx[n], Xl[n], Xr[n], X[n], l,r))

	return Tx[:n+1],X[:n+1], Xl[:n+1], Xr[:n+1]
	
	
	
def create_hist(exp1, exp2, exp3, exp4, tau1, tau2, tau3, tau4, T, c, N, init =0):
	"""
	return an N array X for N realizations of the animal, variances and means
	"""
	
	X = np.zeros(N)
	Y = np.zeros(N)
	if init==2:
		init = random.randint(0,1)

	#create random seeds to get integer seeds
	for i in np.arange(N):	
		S0 = random.randint(1,41234)
		S1 = random.randint(1,9934)
		S2 = random.randint(1,94322)
		S3 = random.randint(1,6193)
		
		t_move1 = get_memory_string(exp1, exp2, S0, S1, tau1, tau2, T, init = init)[1]
		
		t_move2 = get_memory_string(exp3, exp4, S2, S3, tau3, tau4, T, init = init)[1]
		
		#print("{:3.4f} and {:3.4f} \n".format(t_move1, t_move2))
		
		# we dont have to follow the whole path of the animal
		# it is enough to find out how long the animal moved in total 
		# and subtract both times for the respective legs
		X[i] =  -c*t_move1 
		Y[i] = c* t_move2
	
	return np.var(X)+np.mean(X)**2, np.var(Y) + np.mean(Y)**2, np.mean(X), np.mean(Y), np.var(X + Y), X
	

if __name__=="__main__":
	main()
	
	
