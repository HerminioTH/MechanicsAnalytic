import numpy as np
import pandas as pd
import sympy as sy
import os

sec = 1.
minute = 60*sec
hour = 60*minute
day = 24*hour
kPa = 1e3
MPa = 1e6
GPa = 1e9

def double_dot(A, B):
	n, m = A.shape
	value = 0.0
	for i in range(n):
		for j in range(m):
			value += A[i,j]*B[j,i]
	return value

def trace(s):
	return s[0,0] + s[1,1] + s[2,2]

def build_stress(s_xx=0, s_yy=0, s_zz=0, s_xy=0, s_xz=0, s_yz=0):
	return np.array([s_xx, s_yy, s_zz, s_xy, s_xz, s_yz])

def voigt2tensor(s):
	return np.array([
			[s[0], s[3], s[4]],
			[s[3], s[1], s[5]],
			[s[4], s[5], s[2]],
		])

class TensorSaver():
	def __init__(self, output_folder, file_name):
		self.file_name = file_name
		self.output_folder = output_folder
		self.tensor_data = {
			"Time": [], "00": [], "01": [], "02": [], "10": [], "11": [], "12": [], "20": [], "21": [], "22": []
		}

	def save_results(self, time_list, eps_list):
		for time, eps in zip(time_list, eps_list):
			self.__record_values(eps, time)
		self.__save()

	def __record_values(self, tensor, t):
		self.tensor_data["Time"].append(t)
		for i in range(3):
			for j in range(3):
				self.tensor_data[f"{i}{j}"].append(tensor[i,j])

	def __save(self):
		if not os.path.exists(self.output_folder):
			os.makedirs(self.output_folder)
		self.df = pd.DataFrame(self.tensor_data)
		self.df.to_excel(os.path.join(self.output_folder, f"{self.file_name}.xlsx"))



class BaseSolution():
	def __init__(self, settings):
		self.__load_time_list(settings)
		self.__build_sigmas(settings)

	def __load_time_list(self, settings):
		self.time_list = np.array(settings["time"])

	def __build_sigmas(self, settings):
		n = len(settings["sigma_xx"])
		sigma_xx = np.array(settings["sigma_xx"]).reshape((1, n))
		sigma_yy = np.array(settings["sigma_yy"]).reshape((1, n))
		sigma_zz = np.array(settings["sigma_zz"]).reshape((1, n))
		sigma_xy = np.array(settings["sigma_xy"]).reshape((1, n))
		sigma_yz = np.array(settings["sigma_yz"]).reshape((1, n))
		sigma_xz = np.array(settings["sigma_xz"]).reshape((1, n))
		self.sigmas = np.concatenate((sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz)).T
		self.sigma_0 = self.sigmas[0]

	def build_stress_increments(self):
		self.d_sigmas = []
		for i in range(1, len(self.sigmas)):
			self.d_sigmas.append(self.sigmas[i] - self.sigmas[i-1])
		self.d_sigmas = np.array(self.d_sigmas)

	def build_matrix(self, E, nu):
		lame = E*nu/((1+nu)*(1-2*nu))
		G = E/(2 +2*nu)
		x = 1
		self.C = np.array([
					[2*G + lame, 	lame, 		lame, 		0.,		0., 	0.],
					[lame, 			2*G + lame, lame, 		0.,		0., 	0.],
					[lame, 			lame, 		2*G + lame, 0., 	0., 	0.],
					[0., 			0., 		0., 		x*G,	0., 	0.],
					[0., 			0., 		0., 		0., 	x*G, 	0.],
					[0., 			0., 		0., 		0., 	0.,		x*G ],
				])
		self.D = np.linalg.inv(self.C)


class Elastic(BaseSolution):
	def __init__(self, settings):
		super().__init__(settings)
		self.__load_properties(settings)
		self.build_stress_increments()
		self.build_matrix(self.E0, self.nu0)

	def __load_properties(self, settings):
		self.E0 = settings["elasticity"]["E"]
		self.nu0 = settings["elasticity"]["nu"]

	def compute_strains(self):
		self.eps = []
		for i in range(len(self.time_list)):
			eps_value = voigt2tensor(np.dot(self.D, self.sigmas[i]))
			self.eps.append(eps_value)
		self.eps = np.array(self.eps)


class ViscoElastic(BaseSolution):
	def __init__(self, settings):
		super().__init__(settings)
		self.__load_properties(settings)
		self.build_stress_increments()
		self.build_matrix(self.E0, self.nu0)

	def __load_properties(self, settings):
		self.E0 = settings["elasticity"]["E"]
		self.nu0 = settings["elasticity"]["nu"]
		self.voigt_E = np.array(settings["viscoelasticity"]["E"])
		self.voigt_eta = np.array(settings["viscoelasticity"]["eta"])

	def A(self, t):
		try:
			value = np.zeros(len(t))
		except:
			value = 0
		for E, eta in zip(self.voigt_E, self.voigt_eta):
			value += (1 - np.exp(-E*t/eta))*self.E0/E
		return value

		# for i in range(len(self.voigt_E)):
		# 	value = (1 - np.exp(-self.voigt_E*t/self.voigt_eta))*self.E0/self.voigt_E
		# return sum(value)

	def compute_strains(self):
		self.eps = []
		shape = self.d_sigmas[0].shape
		for i in range(0, len(self.time_list)):
			values = np.zeros(shape)
			A_list = self.A(self.time_list[i] - self.time_list[:i])
			A_list = A_list.reshape((len(A_list),1))
			values = A_list*self.d_sigmas[0:i]
			soma = np.sum(values, axis=0)
			eps_value = self.A(self.time_list[i])*voigt2tensor(np.dot(self.D, self.sigma_0))
			eps_value += voigt2tensor(np.dot(self.D, soma))
			self.eps.append(eps_value)
		self.eps = np.array(self.eps)

	# def compute_strains(self):
	# 	self.eps = []
	# 	shape = self.d_sigmas[0].shape
	# 	for i in range(len(self.time_list)):
	# 		soma = np.zeros(shape)
	# 		for j in range(i-1):
	# 			soma += self.A(self.time_list[i] - self.time_list[j])*self.d_sigmas[j+1]
	# 		eps_value = self.A(self.time_list[i])*voigt2tensor(np.dot(self.D, self.sigma_0))
	# 		eps_value += voigt2tensor(np.dot(self.D, soma))
	# 		self.eps.append(eps_value)
	# 	self.eps = np.array(self.eps)

	# def compute_strains(self):
	# 	self.eps = []
	# 	for i in range(len(self.time_list)):
	# 		eps_value = self.A(self.time_list[i])*voigt2tensor(np.dot(self.D, self.sigma_0))
	# 		for j in range(i-1):
	# 			eps_value += self.A(self.time_list[i] - self.time_list[j])*voigt2tensor(np.dot(self.D, self.d_sigmas[j+1]))
	# 		self.eps.append(eps_value)
	# 	self.eps = np.array(self.eps)

class Creep(BaseSolution):
	def __init__(self, settings):
		super().__init__(settings)
		self.__load_properties(settings)

		self.eps_cr = np.zeros((3,3))
		self.eps_cr_old = np.zeros((3,3))
		self.eps_cr_rate = np.zeros((3,3))

	def __load_properties(self, settings):
		self.A = settings["creep"]["A"]
		self.n = settings["creep"]["n"]
		self.T = settings["creep"]["T"]
		self.R = 8.32		# Universal gas constant
		self.Q = 51600  	# Creep activation energy, [J/mol]
		self.B = self.A*np.exp(-self.Q/(self.R*self.T))

	def update_internal_variables(self):
		self.eps_cr_old = self.eps_cr
		self.eps_cr_rate_old = self.eps_cr_rate

	def compute_eps_cr_rate(self, sigma):
		stress = voigt2tensor(sigma)
		s = stress - (1./3)*trace(stress)*np.eye(3)
		von_Mises = np.sqrt((3/2.)*double_dot(s, s))
		self.eps_cr_rate = self.B*(von_Mises**(self.n-1))*s

	def compute_eps_cr(self, i):
		t = self.time_list[i]
		t_old = self.time_list[i-1]
		dt = t - t_old
		self.compute_eps_cr_rate(self.sigmas[i])
		self.eps_cr = self.eps_cr_old + self.eps_cr_rate*dt
		self.eps_cr_old = self.eps_cr

	def compute_strains(self):
		self.eps = [self.eps_cr]
		for i in range(1, len(self.time_list)):
			self.compute_eps_cr(i)
			self.eps.append(self.eps_cr)
		self.eps = np.array(self.eps)





class ViscoPlastic_Desai(BaseSolution):
	def __init__(self, settings):
		super().__init__(settings)
		self.__load_properties(settings)
		self.__initialize_variables()
		self.__initialize_potential_function()
		self.qsi = 0.0

	def compute_strains(self):
		self.eps = [np.zeros((3,3))]
		for i in range(1, len(self.time_list)):
			dt = self.time_list[i] - self.time_list[i-1]
			# dt /= day
			stress_MPa = self.sigmas[i,:].copy()/MPa
			self.compute_yield_function(stress_MPa)
			if self.Fvp <= 0:
				self.eps.append(self.eps[-1])
				self.alphas.append(self.alpha)
				self.alpha_qs.append(self.alpha_q)
			else:
				tol = 1e-6
				error = 2*tol
				maxiter = 50
				alpha_last = self.alpha
				ite = 1
				while error > tol and ite < maxiter:
					strain_rate = self.__compute_strain_rate(stress_MPa)

					increment = double_dot(strain_rate, strain_rate)**0.5*dt
					self.qsi = self.qsi_old + increment

					self.__update_alpha()
					self.__update_alpha_q()

					error = abs(self.alpha - alpha_last)
					alpha_last = self.alpha
					self.compute_yield_function(stress_MPa)

					ite += 1
					if ite >= maxiter:
						print(f"Maximum number of iterations ({maxiter}) reached.")

				self.qsi_old = self.qsi
				self.eps.append(self.eps[-1] + strain_rate*dt)
				self.alphas.append(self.alpha)
				self.alpha_qs.append(self.alpha_q)

		self.eps = np.array(self.eps)
		self.alphas = np.array(self.alphas)
		self.alpha_qs = np.array(self.alpha_qs)

	def __update_alpha(self):
		self.alpha = self.a_1 / (self.qsi**self.eta)

	def __update_alpha_q(self):
		self.alpha_q = self.alpha

	def __compute_strain_rate(self, stress_MPa):
		self.__compute_F0(stress_MPa)
		n_flow = self.evaluate_flow_direction(stress_MPa, self.alpha_q)
		lmbda = self.mu_1*(self.Fvp/self.F0)**self.N_1
		strain_rate = lmbda*n_flow
		return strain_rate

	def __compute_F0(self, stress_MPa):
		# I1, I2, I3 = self.__compute_stress_invariants(*stress_MPa)
		# J1, J2, J3 = self.__compute_deviatoric_invariants(I1, I2, I3)
		# cos3theta = -(J3*np.sqrt(27))/(2*J2**1.5)
		# print("cos3theta: ", cos3theta)
		# print("alpha: ", self.alpha)
		# print("I1: ", I1)
		# print("gamma: ", self.gamma)
		# print("beta_1: ", self.beta_1)
		# print("beta: ", self.beta)
		# print(np.exp(self.beta_1*self.sigma_t))
		# print(np.exp(self.beta*cos3theta))
		# print(np.exp(self.beta_1*self.sigma_t) - self.beta*cos3theta)
		# F1 = (self.gamma*self.sigma_t**2 - self.alpha*self.sigma_t**self.n)
		# F2 = (np.exp(self.beta_1*self.sigma_t) - self.beta*cos3theta)**self.m
		# self.F0 = F1*F2
		# print("F0: ", self.F0)
		# print("F1: ", F1)
		# print("F2: ", F2)
		pass

	def __load_properties(self, settings):
		self.mu_1 = settings["viscoplastic"]["mu_1"]
		self.N_1 = settings["viscoplastic"]["N_1"]
		self.n = settings["viscoplastic"]["n"]
		self.a_1 = settings["viscoplastic"]["a_1"]
		self.eta = settings["viscoplastic"]["eta"]
		self.beta_1 = settings["viscoplastic"]["beta_1"]
		self.beta = settings["viscoplastic"]["beta"]
		self.m = settings["viscoplastic"]["m"]
		self.gamma = settings["viscoplastic"]["gamma"]
		self.k_v = settings["viscoplastic"]["k_v"]
		self.sigma_t = settings["viscoplastic"]["sigma_t"]
		self.alpha_0 = settings["viscoplastic"]["alpha_0"]
		self.F0 = settings["viscoplastic"]["F_0"]

	def __initialize_variables(self):
		self.alpha = self.alpha_0
		self.alpha_q = self.alpha_0
		self.alphas = [self.alpha]
		self.alpha_qs = [self.alpha_q]
		self.qsi_old = (self.a_1/self.alpha)**(1/self.eta)

	def __compute_stress_invariants(self, s_xx, s_yy, s_zz, s_xy, s_xz, s_yz):
		I1 = (s_xx + s_yy + s_zz + 0*self.sigma_t)
		I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
		I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
		return I1, I2, I3

	def __compute_deviatoric_invariants(self, I1, I2, I3):
		J1 = np.zeros(I1.size) if type(I1) == np.ndarray else 0
		J2 = (1/3)*I1**2 - I2
		J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
		return J1, J2, J3

	def __compute_Sr(self, J2, J3):
		return -(J3*np.sqrt(27))/(2*J2**1.5)
		# return (J3**(1/3))/(J2**0.5)

	def compute_yield_function(self, stress_MPa):
		I1, I2, I3 = self.__compute_stress_invariants(*stress_MPa)
		J1, J2, J3 = self.__compute_deviatoric_invariants(I1, I2, I3)
		if J2 == 0.0:
			self.Fvp = -100
		else:
			Sr = self.__compute_Sr(J2, J3)
			I1_star = I1 + self.sigma_t
			F1 = (-self.alpha*I1_star**self.n + self.gamma*I1_star**2)
			F2 = (np.exp(self.beta_1*I1_star) - self.beta*Sr)**self.m
			self.Fvp = J2 - F1*F2

	def __initialize_potential_function(self):
		# Stress components
		self.s_xx = sy.Symbol("s_xx")
		self.s_yy = sy.Symbol("s_yy")
		self.s_zz = sy.Symbol("s_zz")
		self.s_xy = sy.Symbol("s_xy")
		self.s_xz = sy.Symbol("s_xz")
		self.s_yz = sy.Symbol("s_yz")
		self.a_q = sy.Symbol("a_q")

		I1, I2, I3 = self.__compute_stress_invariants(self.s_xx, self.s_yy, self.s_zz, self.s_xy, self.s_xz, self.s_yz)
		J1, J2, J3 = self.__compute_deviatoric_invariants(I1, I2, I3)

		# Compute Lode's angle
		Sr = self.__compute_Sr(J2, J3)

		# Potential function
		F1 = (-self.a_q*I1**self.n + self.gamma*I1**2)
		F2 = (sy.exp(self.beta_1*I1) - self.beta*Sr)**self.m
		self.Qvp = J2 - F1*F2

		variables = (self.s_xx, self.s_yy, self.s_zz, self.s_xy, self.s_xz, self.s_yz, self.a_q)
		self.dQdSxx = sy.lambdify(variables, sy.diff(self.Qvp, self.s_xx), "numpy")
		self.dQdSyy = sy.lambdify(variables, sy.diff(self.Qvp, self.s_yy), "numpy")
		self.dQdSzz = sy.lambdify(variables, sy.diff(self.Qvp, self.s_zz), "numpy")
		self.dQdSxy = sy.lambdify(variables, sy.diff(self.Qvp, self.s_xy), "numpy")
		self.dQdSxz = sy.lambdify(variables, sy.diff(self.Qvp, self.s_xz), "numpy")
		self.dQdSyz = sy.lambdify(variables, sy.diff(self.Qvp, self.s_yz), "numpy")

	def evaluate_flow_direction(self, stress, alpha_q):
		# Analytical derivatives
		dQdS = np.zeros((3,3))
		dQdS[0,0] = self.dQdSxx(*stress, alpha_q)
		dQdS[1,1] = self.dQdSyy(*stress, alpha_q)
		dQdS[2,2] = self.dQdSzz(*stress, alpha_q)
		dQdS[0,1] = dQdS[1,0] = self.dQdSxy(*stress, alpha_q)
		dQdS[0,2] = dQdS[2,0] = self.dQdSxz(*stress, alpha_q)
		dQdS[1,2] = dQdS[2,1] = self.dQdSyz(*stress, alpha_q)
		return dQdS

	# def evaluate_flow_direction(self, stress, alpha_q):
	# 	# Finite difference derivatives
	# 	dSigma = 0.01
	# 	dQdS = np.zeros(6)
	# 	self.compute_yield_function(stress)
	# 	A = self.Fvp
	# 	for i in range(6):
	# 		s = stress.copy()
	# 		s[i] += dSigma
	# 		self.compute_yield_function(s)
	# 		B = self.Fvp
	# 		dQdS[i] = (B - A)/dSigma
	# 	return voigt2tensor(dQdS)




class ViscoPlastic_VonMises(BaseSolution):
	def __init__(self, settings):
		super().__init__(settings)
		self.__load_properties(settings)

	def __load_properties(self, settings):
		self.yield_stress = settings["vonmises"]["yield_stress"]
		self.eta = settings["vonmises"]["eta"]

	def compute_yield_function(self, sigma):
		stress = voigt2tensor(sigma)
		s = stress - (1./3)*trace(stress)*np.eye(3)
		von_Mises = np.sqrt((3/2.)*double_dot(s, s))
		f = von_Mises - self.yield_stress
		if f >= 0:
			self.yield_stress = von_Mises
		return f

	def compute_f_derivatives(self, sigma):
		dSigma = 0.001
		dFdS = np.zeros(6)
		for i in range(6):
			s = sigma.copy()
			s[i] += dSigma
			dFdS[i] = (self.compute_yield_function(sigma) - self.compute_yield_function(s))/dSigma
		return dFdS

	def ramp(self, f):
		return (f + abs(f))/2.

	def compute_strains(self):
		self.eps = [np.zeros((3,3))]
		for i in range(1, len(self.time_list)):
			f = self.compute_yield_function(self.sigmas[i])
			dFdS = self.compute_f_derivatives(self.sigmas[i])
			gamma = self.ramp(f)/self.eta
			self.eps.append(self.eps[-1] - gamma*voigt2tensor(dFdS))
			# if gamma > 0:
			# 	print(gamma*voigt2tensor(dFdS))
		self.eps = np.array(self.eps)


