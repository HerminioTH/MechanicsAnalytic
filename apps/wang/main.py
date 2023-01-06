import numpy as np
import pandas as pd
import json
import os
import sys
import time
sys.path.append(os.path.join("..", "..", "libs"))
from AnalyticSolutions import *

sec = 1.
minute = 60*sec
hour = 60*minute
day = 24*hour
kPa = 1e3
MPa = 1e6
GPa = 1e9

def read_json(file_name):
	with open(file_name, "r") as j_file:
		data = json.load(j_file)
	return data

def save_json(data, file_name):
	with open(file_name, "w") as f:
	    json.dump(data, f, indent=4)

def write_stresses():
	# Read settings
	settings = read_json("settings.json")

	# Define number of time steps
	n_steps = 200

	# Define time levels
	# t_list = np.linspace(0, 900*hour, n_steps)
	t_list = np.linspace(0, 3400*hour, n_steps)
	settings["time"] = list(t_list)

	# Define stress tensors
	settings["sigma_xx"] = list(np.repeat(0.0, n_steps))
	settings["sigma_yy"] = list(np.repeat(0.0, n_steps))
	settings["sigma_xy"] = list(np.repeat(0.0, n_steps))
	settings["sigma_yz"] = list(np.repeat(0.0, n_steps))
	settings["sigma_xz"] = list(np.repeat(0.0, n_steps))

	t1 = 911*hour
	t2 = 1637.1*hour
	t3 = 2700*hour
	s_zz = []
	for t in t_list:
		if t <= t1:
			s_zz.append(12*MPa)
		elif t1 < t and t <= t2:
			s_zz.append(14*MPa)
		elif t2 < t and t <= t3:
			s_zz.append(16*MPa)
		else:
			s_zz.append(6*MPa)
	settings["sigma_zz"] = s_zz

	# Dump to file
	save_json(settings, "settings.json")

def main():
	# Deine output folder
	output_folder = os.path.join("output")

	# Write stresses, if necessary
	# write_stresses()

	# Read settings
	settings = read_json("settings.json")

	# Initialize models
	model_e = Elastic(settings)
	model_cr = Creep(settings)
	model_ve = ViscoElastic(settings)
	model_vp = ViscoPlastic_Desai(settings)

	# Compute strains
	begin = time.time()
	model_e.compute_strains()
	end = time.time()
	print(f"time eps_e: {end - begin}")

	begin = time.time()
	model_cr.compute_strains()
	end = time.time()
	print(f"time eps_cr: {end - begin}")

	begin = time.time()
	model_ve.compute_strains()
	end = time.time()
	print(f"time eps_ve: {end - begin}")

	begin = time.time()
	model_vp.compute_strains()
	end = time.time()
	print(f"time eps_vp: {end - begin}")

	# Compute total strains
	eps_tot = model_e.eps.copy()
	eps_tot += model_cr.eps.copy()
	eps_tot += model_ve.eps.copy()
	eps_tot += model_vp.eps.copy()

	# eps_tot = model_e.eps + model_cr.eps + model_ve.eps + model_vp.eps

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
	saver_eps_tot.save_results(model_cr.time_list, eps_tot)

	# Save stresses
	stresses = [voigt2tensor(model_e.sigmas[i]) for i in range(len(model_e.time_list))]
	saver_stresses = TensorSaver(output_folder, "stresses")
	saver_stresses.save_results(model_cr.time_list, stresses)

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