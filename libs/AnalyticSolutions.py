import numpy as np
import pandas as pd
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
		self.time_list = settings["time"]

	def __build_sigmas(self, settings):
		n = len(settings["sigma_xx"])
		sigma_xx = np.array(settings["sigma_xx"]).reshape((1, n))
		sigma_yy = np.array(settings["sigma_yy"]).reshape((1, n))
		sigma_zz = np.array(settings["sigma_zz"]).reshape((1, n))
		sigma_xy = np.array(settings["sigma_xy"]).reshape((1, n))
		sigma_yz = np.array(settings["sigma_yz"]).reshape((1, n))
		sigma_xz = np.array(settings["sigma_xz"]).reshape((1, n))
		self.sigmas = np.concatenate((sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_yz, sigma_xz)).T
		self.sigma_0 = self.sigmas[0]

	def build_stress_increments(self):
		self.d_sigmas = []
		for i in range(1, len(self.sigmas)):
			self.d_sigmas.append(self.sigmas[i] - self.sigmas[i-1])

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
		self.voigt_E = settings["viscoelasticity"]["E"]
		self.voigt_eta = settings["viscoelasticity"]["eta"]

	def A(self, t):
		value = 0.0
		for E, eta in zip(self.voigt_E, self.voigt_eta):
			r = E/eta
			value += (1 - np.exp(-r*t))*self.E0/E
		return value

	def compute_strains(self):
		self.eps = []
		for i in range(len(self.time_list)):
			eps_value = self.A(self.time_list[i])*voigt2tensor(np.dot(self.D, self.sigma_0))
			for j in range(i-1):
				eps_value += self.A(self.time_list[i] - self.time_list[j])*voigt2tensor(np.dot(self.D, self.d_sigmas[j+1]))
			self.eps.append(eps_value)
		self.eps = np.array(self.eps)

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

class ViscoPlastic(BaseSolution):
	def __init__(self, settings):
		super().__init__(settings)
		self.__load_properties(settings)
		self.__compute_invariants()
		self.__compute_lode_angle()

	def __load_properties(self, settings):
		self.mu_1 = settings["viscoplastic"]["mu_1"]
		self.N_1 = settings["viscoplastic"]["N_1"]
		self.n = settings["viscoplastic"]["n"]
		self.alpha_1 = settings["viscoplastic"]["alpha_1"]
		self.eta_1 = settings["viscoplastic"]["eta_1"]
		self.beta_1 = settings["viscoplastic"]["beta_1"]
		self.beta = settings["viscoplastic"]["beta"]
		self.m_v = settings["viscoplastic"]["m_v"]
		self.gamma = settings["viscoplastic"]["gamma"]
		self.k_v = settings["viscoplastic"]["k_v"]
		self.sigma_t = settings["viscoplastic"]["sigma_t"]
		self.F0 = 1*MPa
		self.alpha_0 = 1*MPa

	def __compute_invariants(self):
		stress = self.sigmas/GPa
		sigma_11 = stress[:,0]
		sigma_22 = stress[:,1]
		sigma_33 = stress[:,2]
		sigma_12 = stress[:,3]
		sigma_23 = stress[:,4]
		sigma_13 = stress[:,5]
		self.I1 = sigma_11 + sigma_22 + sigma_33
		self.I2 = sigma_11*sigma_22 + sigma_22*sigma_33 + sigma_11*sigma_33
		self.I2 += - sigma_12**2 - sigma_23**2 - sigma_13**2
		self.I3 = sigma_11*sigma_22*sigma_33 + 2*sigma_12*sigma_23*sigma_13
		self.I3 += - sigma_33*sigma_12**2 - sigma_11*sigma_23**2 - sigma_22*sigma_13**2
		self.J1 = np.zeros(self.sigmas[:,0].size)
		self.J2 = (1/3)*self.I1**2 - self.I2
		self.J3 = (2/27)*self.I1**3 - (1/3)*self.I1*self.I2 + self.I3

	def __compute_lode_angle(self):
		self.angle = np.zeros(self.J3.size)
		for i in range(self.J3.size):
			self.angle[i] = (1/3)*np.arccos((self.J3[i]**3)/(2*self.J2[i]**1.5))
			print(i, self.J3[i], self.J2[i], (self.J3[i]**3)/(2*self.J2[i]**1.5), self.angle[i])

	# def compute_yield_function(self):

class ViscoPlastic_VonMises(BaseSolution):
	def __init__(self, settings):
		super().__init__(settings)
		self.__load_properties(settings)

	def __load_properties(self, settings):
		self.yield_stress = settings["vonmises"]["yield_stress"]
		self.eta = settings["vonmises"]["eta"]

	def __compute_invariants(self, stress):
		sigma_11 = stress[0]
		sigma_22 = stress[1]
		sigma_33 = stress[2]
		sigma_12 = stress[3]
		sigma_23 = stress[4]
		sigma_13 = stress[5]

		self.I1 = sigma_11 + sigma_22 + sigma_33
		self.I2 = sigma_11*sigma_22 + sigma_22*sigma_33 + sigma_11*sigma_33
		self.I2 += - sigma_12**2 - sigma_23**2 - sigma_13**2
		self.I3 = sigma_11*sigma_22*sigma_33 + 2*sigma_12*sigma_23*sigma_13
		self.I3 += - sigma_33*sigma_12**2 - sigma_11*sigma_23**2 - sigma_22*sigma_13**2
		self.J1 = np.zeros(self.sigmas[:,0].size)
		self.J2 = (1/3)*self.I1**2 - self.I2
		self.J3 = (2/27)*self.I1**3 - (1/3)*self.I1*self.I2 + self.I3

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
			if gamma > 0:
				print(gamma*voigt2tensor(dFdS))
		self.eps = np.array(self.eps)


