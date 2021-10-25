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


def func(x,a, b):
	return a*x**b 	

def incgamma(a,b):
	"""
	returns incomplete (unregularized) gamma function even with negative values a > -1
	"""
	if a < 0:
		return (incgamma(a+1,b) - b**(a) *np.exp(-b) )/a
	else:
		return gammaincc(a,b)*gamma(a)	# while a = 0 exists, the function 
							# has problems, because gamma(0) = inf
		
def round_low_int(z):
	"""
	returns a number rounded to lowest integer
	"""
	if z < np.round(z):
		x = np.round(z) - 1
	else:
		x = np.round(z)
	
	return x
		
def gav_steh_coeff1( n ):
	"""
	returns Gaver Stehfest coefficients for inverse laplace
	ATTENTION: coeff starts with 0, ends with 2n (index)
	"""
	A =np.zeros(2*n+1) # A0 A1 ... A2n
	for k in np.arange(1,2*n+1): # 0.element ist unbesetzt, np.arange geht nur bis 2n
		for j in np.arange(round_low_int((k+1)/2), min(k,n)+1 ):
			A[k] = A[k] + (-1)**(n+k)/np.math.factorial(n)  * j**(n+1) * binom(n,j) * binom(2*j, j) * binom(j, k-j)
	return A

def gav_steh(t,  num, *args):
	"""
	returns the inverse laplace approximation
	"""
	n = 6
	N = len(t) #for time array
	
	A = np.zeros(2*n+1) # gaver stehffest coefficients
	fun = np.zeros(N) # array for function
	
	
	A =gav_steh_coeff1(n) # fill the array
	
	#******************************************************************************************
	# Choose the function

	if num ==0:
		print("Caclulate for 1/x**2. Expect a linear as inverse transformed\n")
		def func(x, *args):
			"""
			[TEST]
			"""
			return 1/x**2
			
	elif num ==1:
		
		
		print("Caclulate the Pareto 2 PDF for alpha = {:.3f} \n".format(args[2] ))
		def func(x, *args):
			"""
			[TEST] returns the Pareto 2 PDF 
			"""
			tau1 = args[0] 
			tau2 = args[1]  
			alpha = args[2] 
			beta = args[3]
			
			return alpha * tau1**alpha * np.exp(x*tau1) * x**alpha *incgamma(-alpha, s*tau1)
		
		def SurPDF(x, exp, tau):
			return (1-Power_law(x, exp, tau))/x

	elif num ==2:
		#print("Caclulate the mean E(x(t)) for alpha = {:.3f} and beta = {:.3f}\n".format(args[2], args[3]))
		
		def Pow_law(x, alpha, tau1):
			"""
			[TEST] returns the Pareto 2 PDF 
			"""
			
			return alpha * tau1**alpha * np.exp(x*tau1) * x**alpha *incgamma(-alpha, s*tau1)

		def func(x, *args):
			"""
			returns exact (analytic) mean for one leg
			"""
			tau1 = args[0] 
			tau2 = args[1]  
			alpha = args[2] 
			beta = args[3] 
			c = args[4]
			return c* (alpha * np.exp(x*tau1)* x**(alpha-2) * tau1**alpha * incgamma(-alpha, x* tau1) *(
			-1 + beta* np.exp(x*tau2)*x**beta *tau2**beta *incgamma(-beta, x*tau2 )   )) / (
			-1 + alpha *beta * np.exp(x*(tau1+tau2)) * x**(alpha +beta) * tau1**alpha * tau2**beta *incgamma(-alpha, x*tau1) * incgamma(-beta, x* tau2))
						
	elif num ==3:
		#print("Caclulate the MSD E(x^2(t)) for alpha = {:.3f} and beta = {:.3f}\n".format(args[2], args[3]))
		
		def Pow_law(x, alpha, tau1):
			"""
			[TEST] returns the Pareto 2 PDF 
			"""
			
			return alpha * tau1**alpha * np.exp(x*tau1) * x**alpha *incgamma(-alpha, s*tau1)

		def func(x, *args):
			"""
			#returns exact (analytic) MSD for one leg
			"""
			tau1 = args[0] 
			tau2 = args[1]  
			alpha = args[2] 
			beta = args[3] 
			c = args[4]
			init = args[5]
			
			a = alpha
			b = beta
			
			if init ==0:
	
				return 2 *c**2* a * np.exp(s*tau1) * (s*tau1)**a * incgamma(-a,s * tau1) * \
			( 1 - b + b * np.exp(s * tau2) * (s * tau2)**b * (-1 + b + s * tau2) *incgamma(-b, s*tau2) + \
			a * b * np.exp(s * tau1) * (s*tau1)**a * incgamma(-a, s*tau1) * \
			(1- np.exp(s*tau2) * (s * tau2)**b * incgamma(-b, s * tau2) * (1+b + s * tau2) + \
			 b * np.exp(2*s * tau2)  * (s * tau2)**(2*b) * incgamma(-b, s*tau2)**2))/ \
			(s**3 * (-1 + a*b * np.exp(s *(tau1 + tau2)) *incgamma(-a, s*tau1) * incgamma(-b, s*tau2)* \
			(s*tau1)**a * (s*tau2)**b )**2)

			elif init == 1:
			# if moving is initial state this is the second moment
				return 2*c**2* (1-b + b * np.exp(s*tau2 ) * (s*tau2)**b * \
				(-1+b + s*tau2)*incgamma(-b, s*tau2) +\
			a * b * np.exp(s * tau1) * (s*tau1)**a * incgamma(-a, s*tau1) * \
			(1 - np.exp(s*tau2) * (s * tau2)**b * incgamma(-b, s * tau2) * (1+b + s * tau2) + \
			 b * np.exp(2*s * tau2)  * (s * tau2)**(2*b) * incgamma(-b, s*tau2)**2))/ \
			(s**3 * (-1 + a*b * np.exp(s *(tau1 + tau2)) *incgamma(-a, s*tau1) * incgamma(-b, s*tau2)* \
			(s*tau1)**a * (s*tau2)**b )**2)
			
	#******************************************************************************************
	# Main Calculation step

	for time in np.arange(N):
		for k in np.arange(1,2*n + 1): # 1 2 ... 2n
			s = k*np.log(2) / t[time]
			ans = func(s, *args)
			#print("t = {:.3f} & A = {:.3f}\n".format(t[time], A[k]))
			fun[time] = fun[time] + A[k] * func(s, *args)

		fun[time] = fun[time]* np.log(2) /t[time]	
	
	return fun
def L1(s, time,*args):
	"""
	varying function for when  beta < alpha, 1			
	"""
	tau1 = args[0] 
	tau2 = args[1]  
	alpha = args[2] 
	beta = args[3] 
	c = args[4]
	
	if time ==2:
		return s**(alpha-beta) * tau1**alpha * gamma(1-alpha) + \
			+ tau2**beta * gamma(1-beta) - \
			s**(1-beta) * (tau2/(1-beta) + tau1/(1-alpha))
		
def L2(s, time,*args):
	"""
	slowly varying function for when alpha < beta, 1			
	"""
	tau1 = args[0] 
	tau2 = args[1]  
	alpha = args[2] 
	beta = args[3] 
	c = args[4]
	
	if time ==2:
		return s**(beta-alpha) * tau2**beta * gamma(1-beta) + \
			+ tau1**alpha * gamma(1-alpha) - \
			s**(1-alpha) * (tau2/(1-beta) + tau1/(1-alpha))
				
def L3(s, time, *args):
	"""
	slowly varying function for when 1 < beta, alpha			
	"""
	tau1 = args[0] 
	tau2 = args[1]  
	alpha = args[2] 
	beta = args[3] 
	c = args[4]
	
	if time ==2:
		return s**(beta-1) * tau2**beta * gamma(1-beta) + \
			s**(alpha-1) * tau1**alpha * gamma(1-alpha) - \
			(tau2/(1-beta) + tau1/(1-alpha))
			
def psi_long(exp,tau,s):
	"""
	returns the asymptotic PDF for psi(s)
	"""
	return 1 + s * tau / (1-exp)  - gamma(1-exp) * tau**exp * s**exp 

def sec_mom_asy(time, t, *args):
	"""
	returns second moment asymptotics
	"""
	tau1 = args[0] 
	tau2 = args[1]  
	alpha = args[2] 
	beta = args[3] 
	c = args[4]
	
	T1 = tau1**alpha * gamma(1-alpha)
	T2 = tau2**beta *gamma(1-beta)
	
	# short term
	if time==0:
		return 2*alpha/tau1*c**2*t**3

	# long range
	elif time ==2:
		a=alpha
		b=beta
		
		if a > b and b > 1:	# region 1
			return  2*c**2 * t**2/L3(1/t,time,*args)**2 * \
			(tau2**2/(1-b)**2 /gamma(3) + \
			t**(1-b)/gamma(4-b) * tau2**b * gamma(1-b) * (1-b) * L3(1/t,time,*args) + \
			t**(2-2*b) * b* gamma(1-b)**2 *tau2**(2*b)/gamma(5-2*b)) 
		elif b> a and a > 1:	# region 2
			return	2*c**2 * t**2/L3(1/t,time,*args)**2 * \
			(tau2**2/(1-b)**2 /gamma(3) + \
			t**(1-b)/gamma(4-b) * tau2**b * gamma(1-b) * (1-b)  * L3(1/t,time,*args) + \
			t**(2-2*b) *b* gamma(1-b)**2 *tau2**(2*b)/gamma(5-2*b)) 

		elif a > 1 and 1 > b:	# region 3
			return  2*c**2 * t**2 /L1(1/t, time, *args)**2 * \
			( gamma(1-b)**2 *tau2**(2*b)/gamma(3) *b+ \
			1/gamma(3) * tau2**b * gamma(1-b) * (1-b) *L1(1/t, time, *args)   + \
			t**(2*b-2)/gamma(1+2*b) *tau2**2/(1-b)**2 )
	
		elif b > 1 and 1 > a:	
			if a + b < 2:	# region 5
				return  2*c**2 * t**(2+a-b) /L2(1/t,time, *args)**2 * psi_long(a, tau1, 1/t)*  \
				( t**(a-b)* gamma(1-b)**2 * b *tau2**(2*b)/gamma(3+2*a-2*b) * psi_long(a, tau1, 1/t)+ \
				1/gamma(3+a-b) * tau2**b * gamma(1-b) * (1-b) *L2(1/t,time, *args) + \
				t**(a+b-2)/gamma(1+2*a) *tau2**2/(1-b)**2* psi_long(a, tau1, 1/t) )
	 
			elif a + b > 2:	# region 4
				return 2 * c**2 * t**(2*a)/L2(1/t, time, *args)**2 * psi_long(a, tau1, 1/t)* \
				( tau2**2/(1-b)**2 / gamma(1+2*a)  * psi_long(a, tau1, 1/t)+ \
				t**(2-a-b)/gamma(3+a-b) * tau2**b * gamma(1-b) * (1-b)  *L2(1/t, time, *args)  +\
				t**(2-2*b)/gamma(3-2*b+2*a) *gamma(1-b)**2 * b* tau2**(2*b) * psi_long(a, tau1, 1/t)) 
				 	
		elif 1 > a and a > b:	# region 6
			return  2*c**2 * t**2 /L1(1/t, time, *args)**2 * \
			( gamma(1-b)**2 *tau2**(2*b)/gamma(3) *b+ \
			1/gamma(3) * tau2**b * gamma(1-b) * (1-b) *L1(1/t, time, *args)   + \
			t**(2*b-2)/gamma(1+2*b) *tau2**2/(1-b)**2 )
				
		elif 1 > b and b > a:	# region 7
			return  2*c**2 * t**(2+a-b) /L2(1/t,time, *args)**2 * psi_long(a, tau1, 1/t)* \
				( t**(a-b)* gamma(1-b)**2 * b *tau2**(2*b)/gamma(3+2*a-2*b) + \
				1/gamma(3+a-b) * tau2**b * gamma(1-b) * (1-b) *L2(1/t,time, *args) + \
				t**(a+b-2)/gamma(1+2*b) *tau2**2/(1-b)**2 )
		
def first_mom_asy(time, t, *args):
	"""
	returns first moment asymptotics
	"""
	tau1 = args[0] 
	tau2 = args[1]  
	alpha = args[2] 
	beta = args[3] 
	c = args[4]
	
	T1 = tau1**alpha * gamma(1-alpha)
	T2 = tau2**beta *gamma(1-beta)
		
	# short term
	if time==0:
		return alpha*c/tau1 * t**2
	# mid range
	elif time==1:
		if alpha > beta:
			return c**2 * t**2 - c**2*t**(2-alpha)*T1/gamma(3-alpha) 
		elif alpha < beta:
			return 2*c**2*(beta - 1) * (t**2/gamma(3) - t**(2-beta+alpha) *T2/T1 /gamma(3-beta+alpha))

	# long range
	elif time ==2:
		return 
	

def main():
	
	tau1 = 2
	tau2 = tau1
	alpha = 0.7
	beta = 0.3

	c = 2.18

	args = (tau1, tau2, alpha, beta,c)
	
	N = 7000	# number of realizations for pushmi pullyu paths
	
	T = np.arange(10, 3000, 400)

	f = np.zeros(len(T))	
	f2 = np.zeros(len(T))
	f =gav_steh(T, 2, *args)
	f2 = gav_steh(T, 3, *args)

	mean_s = np.zeros(len(T))
	MSD_s = np.zeros(len(T))
	y_err = np.zeros(len(T))
	asymp2mom = sec_mom_asy(2, T, *args)
	asymp1mom = first_mom_asy(1,T, *args)
	
	print("T		mean_s		MSD_s \n")
	print("-------------------------------------\n")
	for i in np.arange(len(T)):
		# imported functions to calculate MSD and mean for one leg
		#mean_s[i] = pushmi_pullyu.create_hist(alpha, beta, alpha, beta, tau1, tau2, tau1, tau2, T[i], c, N)[3] 
		MSD_s[i] = pushmi_pullyu.create_hist(alpha, beta, alpha, beta, tau1, tau2, tau1, tau2, T[i], c, N)[1]
		#f2[i] = f2[i] - f[i]**2
		#y_err[i] = pushmi_pullyu.create_hist(alpha, beta, alpha, beta, tau1, tau2, tau1, tau2, T[i], c, N)[6]**0.5
		print("{:3.2f}	{:3.2f}	 \n".format(T[i], f2[i] ))
	
	# fitting
	popt1, popv1 = curve_fit(func, T, f, bounds=(0,[100,3]))
	popt2, popv2 = curve_fit(func, T, f2 , bounds=(0,[100,3]))

	print("func=\n")
	print( func(T, *popt1))
	"""
	plt.figure(1)
	plt.title("Mean")
	plt.plot(T, f, label = 'gaver - stehfest', marker='+')
	plt.plot(T, mean_s, label = 'simulated by model', marker='.')
	plt.grid('True')
	plt.legend(loc=2)
	"""
	plt.figure(2)
	plt.title("MSD")
	plt.loglog(T, f2, label = 'gaver - stehfest', marker='+')
	plt.loglog(T, asymp2mom, label='asymptotics', marker ='o')
	plt.loglog(T, MSD_s, label='simulation model', marker = '*')
	plt.grid('True')
	plt.legend(loc=2)
	

	"""
	plt.figure(3)
	plt.title(" fit function to get parameter")
	plt.plot(T, func(T, *popt1), 'g--', label = 'fit: a = %5.3f, b =%5.3f' %tuple(popt1))
	plt.plot(T, f, label = 'mean by simulation')
	plt.grid('True')
	plt.legend(loc=2)
	
	plt.figure(4)
	plt.title(" fit function to get parameter")
	plt.plot(T, func(T, *popt2), 'g--', label = 'fit: a = %5.3f, b =%5.3f' %tuple(popt2))
	plt.plot(T, f2, label = 'MSD by simulation')
	plt.grid('True')
	plt.legend(loc=2)
	"""
	plt.show()

if __name__=='__main__':
	main()
