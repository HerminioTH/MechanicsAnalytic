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
	data = pd.read_excel(os.path.join("output", "test", "eps_tot.xlsx"))
	print(data)

	time = data["Time"].values/day
	eps_axial = data["22"].values*100
	eps_radial = data["00"].values*100

	fig, ax1 = plt.subplots(1, 1, figsize=(5, 4))
	fig.subplots_adjust(top=0.935, bottom=0.125, left=0.135, right=0.985, hspace=0.2, wspace=0.2)

	ax1.plot(time, eps_axial, ".", color="steelblue", label=r"$\varepsilon_{axial}$")
	ax1.plot(time, eps_radial, ".", color="lightcoral", label=r"$\varepsilon_{radial}$")
	ax1.set_xlabel("Time [day]", size=12, fontname="serif")
	ax1.set_ylabel("Strain [%]", size=12, fontname="serif")
	ax1.grid(True)
	ax1.legend(loc=0, shadow=True, fancybox=True)

	apply_dark_theme(fig, [ax1], transparent=False)

	plt.show()

if __name__ == '__main__':
	main()