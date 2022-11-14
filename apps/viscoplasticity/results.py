import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import os

sec = 1.
minute = 60*sec
hour = 60*minute
day = 24*hour
kPa = 1e3
MPa = 1e6
GPa = 1e9

def apply_grey_theme(fig, axes, transparent=True):
	# fig.patch.set_facecolor("#212121ff")
	# if transparent:
	# 	fig.patch.set_alpha(0.0)
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
	# fig.patch.set_facecolor("#212121ff")
	# if transparent:
	# 	fig.patch.set_alpha(0.0)
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

def main():
	data_tot = pd.read_excel(os.path.join("output", "test_2", "eps_tot.xlsx"))
	data_e = pd.read_excel(os.path.join("output", "test_2", "eps_e.xlsx"))
	data_ve = pd.read_excel(os.path.join("output", "test_2", "eps_ve.xlsx"))
	data_cr = pd.read_excel(os.path.join("output", "test_2", "eps_cr.xlsx"))
	data_vp = pd.read_excel(os.path.join("output", "test_2", "eps_vp.xlsx"))

	time = data_tot["Time"].values/day
	eps_tot_a = data_tot["22"].values*100
	eps_tot_r = data_tot["00"].values*100
	eps_e_a = data_e["22"].values*100
	eps_e_r = data_e["00"].values*100
	eps_ve_a = data_ve["22"].values*100
	eps_ve_r = data_ve["00"].values*100
	eps_vp_a = data_vp["22"].values*100
	eps_vp_r = data_vp["00"].values*100
	eps_cr_a = data_cr["22"].values*100
	eps_cr_r = data_cr["00"].values*100

	fig, ax1 = plt.subplots(1, 1, figsize=(5, 4))
	fig.subplots_adjust(top=0.935, bottom=0.125, left=0.135, right=0.985, hspace=0.2, wspace=0.2)

	ax1.plot(time, eps_tot_a, "-", color="steelblue", label=r"$\varepsilon_{axial}^{tot}$")
	ax1.plot(time, eps_tot_r, "--", color="steelblue", label=r"$\varepsilon_{radial}^{tot}$")
	ax1.plot(time, eps_vp_a, "-", color="gold", label=r"$\varepsilon_{axial}^{vp}$")
	ax1.plot(time, eps_vp_r, "--", color="gold", label=r"$\varepsilon_{radial}^{vp}$")
	# ax1.plot(time, eps_ve_a, "-", color="lightcoral", label=r"$\varepsilon_{axial}^{ve}$")
	# ax1.plot(time, eps_ve_r, "--", color="lightcoral", label=r"$\varepsilon_{radial}^{ve}$")
	# ax1.plot(time, eps_cr_a, "-", color="limegreen", label=r"$\varepsilon_{axial}^{cr}$")
	# ax1.plot(time, eps_cr_r, "--", color="limegreen", label=r"$\varepsilon_{radial}^{cr}$")
	ax1.set_xlabel("Time [day]", size=12, fontname="serif")
	ax1.set_ylabel("Strain [%]", size=12, fontname="serif")
	ax1.grid(True)
	ax1.legend(loc=0, shadow=True, fancybox=True, ncol=2)

	apply_dark_theme(fig, [ax1], transparent=False)

	plt.show()

if __name__ == '__main__':
	main()