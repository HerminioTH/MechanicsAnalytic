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
	sigma_r = 10.0*MPa
	sigma_z = [sigma_i]
	t_mid_1 = t_f/3
	t_mid_2 = 2*t_f/3

	dSigma = 10*MPa
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

	# sigma_r = list(0.5*np.array(sigma_z))

	settings["sigma_zz"] = sigma_z
	settings["sigma_xx"] = list(np.repeat(sigma_r, n_steps))
	settings["sigma_yy"] = list(np.repeat(sigma_r, n_steps))
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
