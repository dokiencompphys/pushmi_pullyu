#!usr/bin/python

"""
This program creates a lot of waiting time histories to test for irregularities. When using the program I encountered a lot of weird drawings due to huge numbers (because of pareto 2 and its fat tails).

"""

import random
import numpy as np
from pushmi_pullyu import get_memory_string 

def main():

	s1= 12345
	s2 = 9876
	N = 8000
	T = 4000
	tau1=2
	tau2=tau1
	
	init=1
	
	alpha= np.random.random_sample(N)+0.01
	beta = np.random.random_sample(N)+0.01

	print("i	alpha		beta		t_move		trest+ tmove		sum all waiting times		observation time")
	for i in range(N):
		s1=random.randint(1,N)
		s2=random.randint(1,N)
		
		t_move =get_memory_string(alpha[i],beta[i],s1, s2, tau1, tau2, T, init=init)[1]
		t1=get_memory_string(alpha[i],beta[i],s1, s2, tau1, tau2, T, init=init)[2]
		t2=get_memory_string(alpha[i],beta[i],s1, s2, tau1, tau2, T, init = init)[4]
		print("{:4.0f}			{:1.2f}			{:1.2f}			{:4.0f}			{:3.2f}		{:3.2f}		{:3.2f}\n".format(i, alpha[i], beta[i],t_move, t1, t2, T))
		
		if np.abs(T-t1) > 10**(-10) or np.abs(T-t2) > 10**(-10):
			print("Test not succesful. Broke up after {:4.0f} runs. alpha = {:1.3f} and beta = {:1.3f}\n".format(i, alpha[i], beta[i]))
			print("T-t1 = {:1e} T-t2= {:1e}\n".format(T-t1, T-t2))
			print(get_memory_string(alpha[i],beta[i],s1, s2, tau1, tau2, T, init = init)[0])
			break
			
	if i==N-1:
		print("Test succesful")

if __name__=='__main__':
	main()
