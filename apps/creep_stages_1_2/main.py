import numpy as np
import pandas as pd
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
	n_steps = 99

	# Define time levels
	settings["time"] = list(np.linspace(0, 150*day, n_steps))

	# Define stress tensors
	settings["sigma_zz"] = list(np.repeat(15*MPa, int(n_steps/3))) + list(np.repeat(5*MPa, int(2*n_steps/3)))
	# settings["sigma_zz"] = list(np.repeat(15*MPa, n_steps))
	settings["sigma_xx"] = list(np.repeat(5*MPa, n_steps))
	settings["sigma_yy"] = list(np.repeat(5*MPa, n_steps))
	settings["sigma_xy"] = list(np.repeat(0.0, n_steps))
	settings["sigma_yz"] = list(np.repeat(0.0, n_steps))
	settings["sigma_xz"] = list(np.repeat(0.0, n_steps))

	# Dump to file
	save_json(settings, "settings.json")

def main():
	# Deine output folder
	output_folder = os.path.join("output", "test")

	# Write stresses, if necessary
	write_stresses()

	# Read settings
	settings = read_json("settings.json")

	# Initialize viscoelastic model
	model_cr = Creep(settings)
	model_ve = ViscoElastic(settings)

	# Compute strains
	model_cr.compute_strains()
	model_ve.compute_strains()

	# Save results
	saver_eps_cr = TensorSaver(output_folder, "eps_cr")
	saver_eps_cr.save_results(model_cr.time_list, model_cr.eps)

	saver_eps_ve = TensorSaver(output_folder, "eps_ve")
	saver_eps_ve.save_results(model_ve.time_list, model_ve.eps)

	saver_eps = TensorSaver(output_folder, "eps_tot")
	saver_eps.save_results(model_ve.time_list, model_ve.eps + model_cr.eps)

	# Save settings
	save_json(settings, os.path.join(output_folder, "settings.json"))




if __name__ == '__main__':
	main()
