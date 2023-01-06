import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import sympy as sy
import os
import json

sec = 1.
minute = 60*sec
hour = 60*minute
day = 24*hour
Pa = 1
kPa = 1e3
MPa = 1e6
GPa = 1e9

def read_json(file_name):
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def apply_grey_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.spines['bottom'].set_color('white')
			ax.spines['top'].set_color('white')
			ax.spines['right'].set_color('white')
			ax.spines['left'].set_color('white')
			ax.tick_params(axis='x', colors='white', which='both')
			ax.tick_params(axis='y', colors='white', which='both')
			ax.yaxis.label.set_color('white')
			ax.xaxis.label.set_color('white')
			ax.title.set_color('white')
			ax.set_facecolor("#2b2b2bff")
			ax.grid(True, color='#414141ff')

def apply_dark_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#37474fff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.spines['bottom'].set_color('white')
			ax.spines['top'].set_color('white')
			ax.spines['right'].set_color('white')
			ax.spines['left'].set_color('white')
			ax.tick_params(axis='x', colors='white', which='both')
			ax.tick_params(axis='y', colors='white', which='both')
			ax.yaxis.label.set_color('white')
			ax.xaxis.label.set_color('white')
			ax.title.set_color('white')
			ax.set_facecolor("#424242ff")
			ax.grid(True, color='#5a5a5aff')

def apply_white_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
		if ax != None:
			ax.grid(True, color='#c4c4c4ff')
			ax.set_axisbelow(True)
			ax.spines['bottom'].set_color('black')
			ax.spines['top'].set_color('black')
			ax.spines['right'].set_color('black')
			ax.spines['left'].set_color('black')
			ax.tick_params(axis='x', colors='black', which='both')
			ax.tick_params(axis='y', colors='black', which='both')
			ax.yaxis.label.set_color('black')
			ax.xaxis.label.set_color('black')
			ax.set_facecolor("#e9e9e9ff")

def load_stresses(settings, unit=MPa):
	s_xx = np.array(settings["sigma_xx"])/unit
	s_yy = np.array(settings["sigma_yy"])/unit
	s_zz = np.array(settings["sigma_zz"])/unit
	s_xy = np.array(settings["sigma_xy"])/unit
	s_yz = np.array(settings["sigma_yz"])/unit
	s_xz = np.array(settings["sigma_xz"])/unit
	return s_xx, s_yy, s_zz, s_xy, s_yz, s_xz

def compute_invariants(settings, s_xx, s_yy, s_zz, s_xy, s_yz, s_xz):
	sigma_t = settings["viscoplastic"]["sigma_t"]
	I1 = (s_xx + s_yy + s_zz + 0*sigma_t)
	I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
	I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
	J1 = np.repeat(0., len(I1))
	J2 = (1/3)*I1**2 - I2
	J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
	Sr = []
	for j2, j3 in zip(J2, J3):
		Sr.append(-(j3*np.sqrt(27))/(2*j2**1.5))
	return I1, np.sqrt(J2), np.array(Sr)

def plot_yield_surface(ax, Sr, alpha, settings):
	gamma = settings["viscoplastic"]["gamma"]
	n = settings["viscoplastic"]["n"]
	m = settings["viscoplastic"]["m"]
	beta = settings["viscoplastic"]["beta"]
	beta_1 = settings["viscoplastic"]["beta_1"]
	sigma_t = settings["viscoplastic"]["sigma_t"]#/MPa

	# Yield surface
	I1 = np.linspace(-sigma_t, 40, 1500)
	I1_star = I1 + sigma_t
	J2 = (-alpha*I1_star**n + gamma*I1_star**2)*(np.exp(beta_1*I1_star) - beta*Sr)**m
	np.clip(J2, a_min=0.0, a_max=None, out=J2)

	# Dilatancy boundary
	J2_dil = (1 - 2/n)*gamma*(I1_star**2)*(np.exp(beta_1*I1_star) - beta*Sr)**m

	# Short-term failure boundary
	alpha = 0
	J2_fail = (-alpha*I1_star**n + gamma*I1_star**2)*(np.exp(beta_1*I1_star) - beta*Sr)**m

	ax.plot(I1, np.sqrt(J2), "-", color="steelblue")
	ax.plot(I1, np.sqrt(J2_dil), "--", color="lightcoral")
	ax.plot(I1, np.sqrt(J2_fail), "--", color="forestgreen")
	# ax.set_ylim(0, 25)

def read_data_exp():
	data = pd.read_pickle(os.path.join("..", "..", "..", "..", "..", "experimental", "experiment_BBG_Pr-5", "experiment_BBG_Pr-5.pkl"))
	data = data.iloc[1817:]
	time = data["Experiment Time (hours)"].values
	eps_axial = data["Axial Strain (dL/Lo)"].values*100
	eps_radial = data["Radial Strain (dR/Ro)"].values*100
	time = time - time[0]
	eps_axial = eps_axial - eps_axial[0]
	eps_radial = eps_radial - eps_radial[0]
	stress_axial = data["Axial Eff. Stress (MPa)"].values - data["Axial Eff. Stress (MPa)"].values[0]
	stress_radial = data["Radial Eff. Stress (MPa)"].values - data["Radial Eff. Stress (MPa)"].values[0]
	return time, eps_axial, eps_radial, stress_axial, stress_radial

def plot_exp(ax):
	data = np.loadtxt(os.path.join("data", "data.csv"), delimiter=",")
	time = data[:,0]
	eps = data[:,1]
	ax.plot(time, eps, ".", color="steelblue", label="Wang")

def plot_surfaces(ax1, settings, alphas, stresses_num):
	# Compute invariants
	I1, sqrt_J2, Sr_list = compute_invariants(settings, *stresses_num)

	for Sr, alpha in zip(Sr_list[0:], alphas[0:]):
		plot_yield_surface(ax1, Sr, alpha, settings)
	ax1.plot(I1, sqrt_J2, ".", color="gold")
	ax1.set_xlabel(r"$I_1$ (MPa)", fontsize=12, fontname="serif")
	ax1.set_ylabel(r"$\sqrt{J_2}$ (MPa)", fontsize=12, fontname="serif")
	ax1.grid(True)

def main():
	# Output folder
	output_folder = os.path.join("output")

	# Read settings
	settings = read_json("settings.json")

	# Read numerical stresses
	stresses_num = load_stresses(settings, MPa)

	# Load alphas for desai failure model
	data_alpha = pd.read_excel(os.path.join(output_folder, "alpha.xlsx"))
	alphas = data_alpha["00"]

	# print(*stresses_num)

	# Load total strain
	data_tot = pd.read_excel(os.path.join(output_folder, "eps_tot.xlsx"))
	data_e = pd.read_excel(os.path.join(output_folder, "eps_e.xlsx"))
	data_ve = pd.read_excel(os.path.join(output_folder, "eps_ve.xlsx"))
	data_cr = pd.read_excel(os.path.join(output_folder, "eps_cr.xlsx"))
	data_vp = pd.read_excel(os.path.join(output_folder, "eps_vp.xlsx"))

	time = data_tot["Time"].values/hour
	eps_axial_tot = data_tot["22"].values*100
	eps_radial_tot = data_tot["00"].values*100
	eps_vol_tot = eps_axial_tot + 2*eps_radial_tot
	eps_axial_e = data_e["22"].values*100
	eps_axial_ve = data_ve["22"].values*100
	eps_axial_cr = data_cr["22"].values*100
	eps_axial_vp = data_vp["22"].values*100
	eps_radial_vp = data_vp["00"].values*100

	fig, axes = plt.subplots(2, 2, figsize=(14, 8))
	fig.subplots_adjust(top=0.935, bottom=0.125, left=0.065, right=0.985, hspace=0.2, wspace=0.275)
	ax1 = axes[0,0]
	ax2 = axes[0,1]
	ax3 = axes[1,0]
	ax4 = axes[1,1]

	plot_surfaces(ax1, settings, alphas, stresses_num)

	plot_exp(ax2)
	ax2.plot(time, eps_axial_tot, "--", color="steelblue", label="Model")
	ax2.set_xlabel("Time [hours]", size=12, fontname="serif")
	ax2.set_ylabel("Strain [%]", size=12, fontname="serif")
	ax2.grid(True)
	ax2.legend(loc=0, shadow=True, fancybox=True)

	ax3.plot(time, alphas, "-", color="gold", linewidth=2.0)
	ax3.set_xlabel("Time [day]", size=12, fontname="serif")
	ax3.set_ylabel("Alpha", size=12, fontname="serif")

	# ax4.plot(time_exp, stress_axial, "-", color="steelblue", label=r"$\sigma_{axial}$", linewidth=2.0)
	# ax4.plot(time_exp, stress_radial, "-", color="lightcoral", label=r"$\sigma_{radial}$", linewidth=2.0)
	# ax4.set_xlabel("Time [day]", size=12, fontname="serif")
	# ax4.set_ylabel("Stress [MPa]", size=12, fontname="serif")
	# ax4.grid(True)
	# ax4.legend(loc=0, shadow=True, fancybox=True)

	# tao_exp = stress_axial - stress_radial
	# tao_num = stresses_num[2] - stresses_num[0]
	# ax5.plot(eps_axial, tao_exp, "-", color="steelblue")
	# ax5.plot(eps_axial_tot, tao_num, "--", color="steelblue")

	# ax6.plot(eps_radial, tao_exp, "-", color="lightcoral")
	# ax6.plot(eps_radial_tot, tao_num, "--", color="lightcoral")


	apply_dark_theme(fig, axes.flatten(), transparent=False)

	plt.show()

if __name__ == '__main__':
	main()