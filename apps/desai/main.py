import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
import os
import sys
sys.path.append(os.path.join("..", "..", "libs"))
from AnalyticSolutions import *

def read_json(file_name):
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def save_json(data, file_name):
	with open(file_name, "w") as f:
	    json.dump(data, f, indent=4)

def apply_white_theme(fig, axes, transparent=True):
	fig.patch.set_facecolor("#212121ff")
	if transparent:
		fig.patch.set_alpha(0.0)
	for ax in axes:
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

def compute_I1_J2_Lode(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz):
	I1 = s_xx + s_yy + s_zz
	I2 = s_xx*s_yy + s_yy*s_zz + s_xx*s_zz - s_xy**2 - s_yz**2 - s_xz**2
	I3 = s_xx*s_yy*s_zz + 2*s_xy*s_yz*s_xz - s_zz*s_xy**2 - s_xx*s_yz**2 - s_yy*s_xz**2
	J1 = 0
	J2 = (1/3)*I1**2 - I2
	J3 = (2/27)*I1**3 - (1/3)*I1*I2 + I3
	lode_angle = []
	for j2, j3 in zip(J2, J3):
		value = sy.re((1/3)*sy.acos(-(j3*27**0.5)/(2*j2**1.5)))
		lode_angle.append(np.degrees(float(value)))
	return I1, np.sqrt(J2), np.array(lode_angle)

def compute_invariants(settings):
	s_xx = np.array(settings["sigma_xx"])/MPa
	s_yy = np.array(settings["sigma_yy"])/MPa
	s_zz = np.array(settings["sigma_zz"])/MPa
	s_xy = np.array(settings["sigma_xy"])/MPa
	s_yz = np.array(settings["sigma_yz"])/MPa
	s_xz = np.array(settings["sigma_xz"])/MPa
	I1, sqrt_J2, lode = compute_I1_J2_Lode(s_xx, s_yy, s_zz, s_xy, s_xz, s_yz)
	return I1, sqrt_J2, lode

def plot_stress_path(ax, I1, sqrt_J2):
	ax.plot(I1, sqrt_J2, ".", color="forestgreen")
	ax.set_xlabel(r"$I_1$ (MPa)", fontsize=12, fontname="serif")
	ax.set_ylabel(r"$\sqrt{J_2}$ (MPa)", fontsize=12, fontname="serif")
	ax.grid(True)

def plot_yield_surface(ax, lode_angle, settings):
	gamma = settings["viscoplastic"]["gamma"]
	alpha = settings["viscoplastic"]["alpha_0"]
	n = settings["viscoplastic"]["n"]
	m = settings["viscoplastic"]["m"]
	beta = settings["viscoplastic"]["beta"]
	beta_1 = settings["viscoplastic"]["beta_1"]
	sigma_t = settings["viscoplastic"]["sigma_t"]/MPa

	# Yield surface
	theta_rad = np.radians(lode_angle)
	I1 = np.linspace(-sigma_t, 125, 500)
	I1_star = I1 + sigma_t
	J2 = (-alpha*I1_star**n + gamma*I1_star**2)*(np.exp(beta_1*I1_star) - beta*np.cos(3*theta_rad))**m
	np.clip(J2, a_min=0.0, a_max=None, out=J2)

	# Dilatancy boundary
	J2_dil = (1 - 2/n)*gamma*(I1_star**2)*(np.exp(beta_1*I1_star) - beta*np.cos(3*theta_rad))**m

	ax.plot(I1, np.sqrt(J2), "-", color="steelblue")
	ax.plot(I1, np.sqrt(J2_dil), "--", color="lightcoral")
	ax.set_ylim(0, 25)

def plot_yield_surfaces(ax, lode_angle, settings):
	gamma = settings["viscoplastic"]["gamma"]
	alpha = settings["viscoplastic"]["alpha_0"]
	n = settings["viscoplastic"]["n"]
	m = settings["viscoplastic"]["m"]
	beta = settings["viscoplastic"]["beta"]
	beta_1 = settings["viscoplastic"]["beta_1"]
	sigma_t = settings["viscoplastic"]["sigma_t"]/MPa

	# Yield surface
	theta_rad = np.radians(lode_angle)
	I1 = np.linspace(-sigma_t, 125, 500)
	I1_star = I1 + sigma_t
	J2 = (-alpha*I1_star**n + gamma*I1_star**2)*(np.exp(beta_1*I1_star) - beta*np.cos(3*theta_rad))**m
	np.clip(J2, a_min=0.0, a_max=None, out=J2)

	# Dilatancy boundary
	J2_dil = (1 - 2/n)*gamma*(I1_star**2)*(np.exp(beta_1*I1_star) - beta*np.cos(3*theta_rad))**m

	ax.plot(I1, np.sqrt(J2), "-", color="steelblue")
	ax.plot(I1, np.sqrt(J2_dil), "--", color="lightcoral")
	ax.set_ylim(0, 25)




def write_stresses():
	j_file = "settings.json"

	# Read settings
	settings = read_json(j_file)

	# Define number of time steps
	n_steps = 500

	# Define time levels
	t_f = 25000*hour
	settings["time"] = list(np.linspace(0, t_f, n_steps))

	# Define stress tensors
	# sigma_i = 10.0*MPa
	# sigma_f = 20*MPa
	# sigma_r = 10.0*MPa
	# sigma_z = []
	# t_mid_1 = 5000*hour
	# t_mid_2 = t_f - t_mid_1

	# for t in settings["time"]:
	# 	if t <= t_mid_1:
	# 		sigma_z.append(sigma_i + (sigma_f - sigma_i)*t/t_mid_1)
	# 		# sigma_z.append(10.0*MPa)
	# 	elif t_mid_1 < t and t < t_mid_2:
	# 		sigma_z.append(sigma_f)
	# 	else:
	# 		sigma_z.append(sigma_f)

	sigma_i = 10.0*MPa
	sigma_f1 = 15*MPa
	sigma_f2 = 20*MPa
	sigma_r = 5.0*MPa
	sigma_z = [sigma_i]
	t_mid_1 = t_f/3
	t_mid_2 = 2*t_f/3

	dSigma = 20*MPa
	dt = t_f/4
	m0 = dSigma/dt
	m1 = 0
	m2 = dSigma/dt
	m3 = 0*MPa/dt
	m4 = -dSigma/dt
	m5 = 0*MPa/dt
	m6 = 2*dSigma/dt
	m7 = 0*MPa/dt
	m8 = 0*MPa/dt

	for i in range(1, len(settings["time"])):
		t = settings["time"][i]
		dt = settings["time"][i] - settings["time"][i-1]
		if t <= t_f/9:
			sigma_z.append(sigma_z[-1] + m0*dt)
		elif 1*t_f/9 < t and t <= 2*t_f/9:
			sigma_z.append(sigma_z[-1] + m1*dt)
		elif 2*t_f/9 < t and t <= 3*t_f/9:
			sigma_z.append(sigma_z[-1] + m2*dt)
		elif 3*t_f/9 < t and t <= 4*t_f/9:
			sigma_z.append(sigma_z[-1] + m3*dt)
		elif 4*t_f/9 < t and t <= 5*t_f/9:
			sigma_z.append(sigma_z[-1] + m4*dt)
		elif 5*t_f/9 < t and t <= 6*t_f/9:
			sigma_z.append(sigma_z[-1] + m5*dt)
		elif 6*t_f/9 < t and t <= 7*t_f/9:
			sigma_z.append(sigma_z[-1] + m6*dt)
		elif 7*t_f/9 < t and t <= 8*t_f/9:
			sigma_z.append(sigma_z[-1] + m7*dt)
		else:
			sigma_z.append(sigma_z[-1] + m8*dt)

	sigma_r = list(0.5*np.array(sigma_z))

	settings["sigma_zz"] = sigma_z
	settings["sigma_xx"] = sigma_r#list(np.repeat(sigma_r, n_steps))
	settings["sigma_yy"] = sigma_r#list(np.repeat(sigma_r, n_steps))
	settings["sigma_xy"] = list(np.repeat(0.0*MPa, n_steps))
	settings["sigma_yz"] = list(np.repeat(0.0*MPa, n_steps))
	settings["sigma_xz"] = list(np.repeat(0.0*MPa, n_steps))

	# Dump to file
	save_json(settings, "settings.json")

def main():
	# Deine output folder
	output_folder = os.path.join("output", "test_0")

	# Write stresses, if necessary
	write_stresses()

	# Read settings
	settings = read_json("settings.json")

	# Compute invariants
	I1, sqrt_J2, lode = compute_invariants(settings)

	# Plot stress path and yield surface
	fig, ax = plt.subplots(1, 1, figsize=(6, 3.6))
	fig.subplots_adjust(top=0.92, bottom=0.135, left=0.105, right=0.985, hspace=0.2, wspace=0.2)

	plot_stress_path(ax, I1, sqrt_J2)
	plot_yield_surface(ax, 60, settings)

	apply_white_theme(fig, [ax], transparent=True)
	# plt.show()


	# Initialize models
	model_e = Elastic(settings)
	model_cr = Creep(settings)
	model_ve = ViscoElastic(settings)
	model_vp = ViscoPlastic_Desai(settings)


	# Compute strains
	model_e.compute_strains()
	model_cr.compute_strains()
	model_ve.compute_strains()
	model_vp.compute_strains()

	# Save results
	saver_eps_e = TensorSaver(output_folder, "eps_e")
	saver_eps_e.save_results(model_e.time_list, model_e.eps)

	saver_eps_ve = TensorSaver(output_folder, "eps_ve")
	saver_eps_ve.save_results(model_ve.time_list, model_ve.eps)

	saver_eps_cr = TensorSaver(output_folder, "eps_cr")
	saver_eps_cr.save_results(model_cr.time_list, model_cr.eps)

	saver_eps_vp = TensorSaver(output_folder, "eps_vp")
	saver_eps_vp.save_results(model_vp.time_list, model_vp.eps)

	saver_eps_tot = TensorSaver(output_folder, "eps_tot")
	eps_tot = model_e.eps
	eps_tot += model_ve.eps
	# eps_tot += model_cr.eps
	eps_tot += model_vp.eps
	saver_eps_tot.save_results(model_cr.time_list, eps_tot)

	# Save alphas for viscoplasticity
	saver_alpha = TensorSaver(output_folder, "alpha")
	alphas = np.zeros(model_e.eps.shape)
	alphas[:,0,0] = model_vp.alphas
	saver_alpha.save_results(model_vp.time_list, alphas)

	saver_alpha_q = TensorSaver(output_folder, "alpha_q")
	alpha_qs = np.zeros(model_e.eps.shape)
	alpha_qs[:,0,0] = model_vp.alpha_qs
	saver_alpha_q.save_results(model_vp.time_list, alpha_qs)

	# Save settings
	save_json(settings, os.path.join(output_folder, "settings.json"))




if __name__ == '__main__':
	main()
