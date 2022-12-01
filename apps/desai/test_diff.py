import sympy as sy
import numpy as np

def main():
	# Stress components
	s_xx = sy.Symbol("s_xx")
	s_yy = sy.Symbol("s_yy")
	s_zz = sy.Symbol("s_zz")
	s_xy = sy.Symbol("s_xy")
	s_xz = sy.Symbol("s_xz")
	s_yz = sy.Symbol("s_yz")

	# Compute invariants
	I1 = s_xx + s_yy + s_zz
	I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
	I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2

	J1 = 0
	J2 = (1/3)*I1**2 - I2
	J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3

	# Compute Lode's angle
	theta = (1/3)*sy.acos(-(J3*np.sqrt(27))/(2*J2**1.5))

	theta_val = sy.re(theta.subs([(s_xx, 0), (s_yy, 0), (s_zz, 2e6), (s_xy, 0), (s_xz, 0), (s_yz, 0)]))

	# print(theta)
	# print()	
	# print(float(theta_val))
	# print(np.degrees(float(theta_val)))

	# Material parameters
	a1 = 0.00005
	eta = 0.7
	beta1 = 4.8e-3
	beta = 0.995
	m = -0.5
	n = 3
	N1 = 3
	gamma = 0.11
	F_0 = 1
	T_a = 5.4

	# Parameter k
	# k = -0.001*s_zz**2 + 0.02*s_zz + 0.06

	# Parameter alpha
	alpha = 0.0
	alpha_q = alpha #+ k*(alpha_0 - alpha)*(1 - qsi_v/qsi)

	Q = J2 - (gamma*I1**2 - alpha_q*I1**n)*(sy.exp(beta1*I1) - beta*sy.cos(3*theta))**m

	# print(Q.evalf(subs={s_xx: 10e6}))
	dQdsxx = sy.diff(Q, s_xx)
	# print(dQdsxx)


	x = sy.Symbol("x")

	f = sy.cos(x**2)
	dfdx = sy.diff(f, x)
	print(dfdx)


if __name__ == '__main__':
	main()