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

def build_stress_increments(sigmas):
	d_sigma_list = []
	for i in range(1, len(sigmas)):
		d_sigma_list.append(sigmas[i] - sigmas[i-1])
	return d_sigma_list

class TensorSaver():
	def __init__(self, name):
		self.name = name
		self.tensor_data = {
			"Time": [],
			"00": [],
			"01": [],
			"02": [],
			"10": [],
			"11": [],
			"12": [],
			"20": [],
			"21": [],
			"22": []
		}

	def record_values(self, tensor, t):
		self.tensor_data["Time"].append(t)
		for i in range(3):
			for j in range(3):
				self.tensor_data[f"{i}{j}"].append(tensor[i,j])

	def save(self, output_folder):
		if not os.path.exists(output_folder):
			os.makedirs(output_folder)
		self.df = pd.DataFrame(self.tensor_data)
		self.df.to_excel(os.path.join(output_folder, f"{self.name}.xlsx"))

class Creep():
	def __init__(self, A, n):
		self.A = A
		self.n = n
		self.T = 298  		# Temperature, [K]
		self.R = 8.32		# Universal gas constant
		self.Q = 51600  	# Creep activation energy, [J/mol]
		self.B = self.A*np.exp(-self.Q/(self.R*self.T))
		self.eps_cr = np.zeros((3,3))
		self.eps_cr_old = np.zeros((3,3))
		self.eps_cr_rate = np.zeros((3,3))
		self.t_old = 0

	def update_internal_variables(self):
		self.eps_cr_old = self.eps_cr
		self.eps_cr_rate_old = self.eps_cr_rate

	def compute_eps_cr_rate(self, sigma):
		stress = voigt2tensor(sigma)
		s = stress - (1./3)*trace(stress)*np.eye(3)
		von_Mises = np.sqrt((3/2.)*double_dot(s, s))
		self.eps_cr_rate = self.B*(von_Mises**(self.n-1))*s

	def compute_eps_cr(self, stress, t):
		self.compute_eps_cr_rate(stress)
		self.eps_cr = self.eps_cr_old + self.eps_cr_rate*(t - self.t_old)
		# self.eps_cr = self.eps_cr_rate*t
		self.t_old = t
		self.eps_cr_old = self.eps_cr

class ViscoElastic():
	def __init__(self, E0, nu0, voigt_E, voigt_eta, sigma_0, d_sigmas, time_list, strain_0=voigt2tensor(build_stress(s_xx=0, s_yy=0, s_zz=0, s_xy=0, s_xz=0, s_yz=0))):
		self.E0 = E0
		self.nu0 = nu0
		self.voigt_E = voigt_E
		self.voigt_eta = voigt_eta
		self.sigma_0 = sigma_0
		self.time_list = time_list
		self.d_sigmas = d_sigmas
		self.strain_0 = strain_0
		self.build_matrix()

	def build_matrix(self):
		lame = self.E0*self.nu0/((1+self.nu0)*(1-2*self.nu0))
		G = self.E0/(2 +2*self.nu0)
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

	def A(self, t):
		value = 0.0
		for E, eta in zip(self.voigt_E, self.voigt_eta):
			r = E/eta
			value += (1 - np.exp(-r*t))*self.E0/E
		return 1 + value

	def compute_strains(self):
		self.eps_ve = []
		for i in range(len(self.time_list)):
			eps = self.A(self.time_list[i])*voigt2tensor(np.dot(self.D, self.sigma_0))
			for j in range(i-1):
				eps += self.A(self.time_list[i] - self.time_list[j])*voigt2tensor(np.dot(self.D, self.d_sigmas[j+1]))
			self.eps_ve.append(self.strain_0 + eps)

class ViscoPlastic():
	def __init__(self, A, n):
		self.A = A
		self.n = n
		self.T = 298  		# Temperature, [K]
		self.R = 8.32		# Universal gas constant
		self.Q = 51600  	# Creep activation energy, [J/mol]
		self.B = self.A*np.exp(-self.Q/(self.R*self.T))
		self.eps_cr = np.zeros((3,3))
		self.eps_cr_old = np.zeros((3,3))
		self.eps_cr_rate = np.zeros((3,3))
		self.t_old = 0

	def update_internal_variables(self):
		self.eps_cr_old = self.eps_cr
		self.eps_cr_rate_old = self.eps_cr_rate

	def compute_eps_vp_rate(self, sigma):
		stress = voigt2tensor(sigma)
		s = stress - (1./3)*trace(stress)*np.eye(3)
		von_Mises = np.sqrt((3/2.)*double_dot(s, s))
		self.eps_cr_rate = self.B*(von_Mises**(self.n-1))*s

	def compute_creep_strain(self, stress, t):
		self.compute_creep_strain_rate(stress)
		self.eps_cr = self.eps_cr_old + self.eps_cr_rate*(t - self.t_old)
		self.t_old = t
		self.eps_cr_old = self.eps_cr