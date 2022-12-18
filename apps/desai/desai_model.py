import sympy as sy
import numpy as np
import matplotlib.pyplot as plt

MPa = 1e6

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

def compute_J2_for_F_equal_0(I1, theta_deg, alpha=0.0, gamma=0.11, n=3, m=-0.5, beta=0.995, beta1=4.8e-9, Ta=5.4):
	theta_rad = np.radians(theta_deg)
	I1_star = I1 + Ta
	J2 = (-alpha*I1_star**n + gamma*I1_star**2)*(np.exp(beta1*I1_star) - beta*np.cos(3*theta_rad))**m
	np.clip(J2, a_min=0.0, a_max=None, out=J2)
	return np.sqrt(J2)

def compute_dilatancy_boundary(I1, theta_deg, gamma=0.11, n=3, m=-0.5, beta=0.995, beta1=4.8e-9, Ta=5.4):
	theta_rad = np.radians(theta_deg)
	I1_star = I1 + Ta
	J2 = (1 - 2/n)*gamma*(I1_star**2)*(np.exp(beta1*I1_star) - beta*np.cos(3*theta_rad))**m
	return np.sqrt(J2)

def plot_failure_curve(ax, alpha, theta_deg, I1_list, gamma=0.11, beta=-0.995, beta1=4.8e-3, Ta=5.4/3):
	J2_sqr = compute_J2_for_F_equal_0(I1_list, alpha=alpha, theta_deg=theta_deg, gamma=gamma, beta=beta, beta1=beta1, Ta=Ta)
	ax.plot(I1_list, J2_sqr, "-", color="steelblue")

def plot_dilatancy_boundary(ax, theta_deg, I1_list, gamma=0.11, beta=-0.995, beta1=4.8e-3, Ta=5.4/3):
	J2_dil = compute_dilatancy_boundary(I1_list, theta_deg=theta_deg, gamma=gamma, beta=beta, beta1=beta1, Ta=Ta)
	ax.plot(I1_list, J2_dil, "--", color="lightcoral")

def main():
	beta = -0.995
	beta1 = 0.00479
	gamma = 0.0945
	Ta = 5.4

	I1_list = np.linspace(-Ta, 125, 500)

	fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))
	fig.subplots_adjust(top=0.935, bottom=0.125, left=0.050, right=0.985, hspace=0.2, wspace=0.2)

	for theta, ax in zip([0, 30, 60], [ax1, ax2, ax3]):
		plot_failure_curve(ax, 0.0, theta, I1_list, gamma, beta, beta1, Ta)
		plot_failure_curve(ax, 0.0002, theta, I1_list, gamma, beta, beta1, Ta)
		plot_failure_curve(ax, 0.0005, theta, I1_list, gamma, beta, beta1, Ta)
		plot_failure_curve(ax, 0.0008, theta, I1_list, gamma, beta, beta1, Ta)
		plot_failure_curve(ax, 0.0016, theta, I1_list, gamma, beta, beta1, Ta)
		plot_dilatancy_boundary(ax, theta, I1_list, gamma, beta, beta1, Ta)
		ax.set_xlabel(r"$I_1$" + " (MPa)")
		ax.set_ylabel(r"$\sqrt{J_2}$" + " (MPa)")
		ax.grid(True)
		ax.set_ylim(0, 45)
		ax.set_title(r"$\theta =$" + str(theta) + r"$^o$")

	plt.show()

def main2():
	beta = -0.995
	beta1 = 0.00479
	gamma = 0.0945
	Ta = 5.4

	I1_list = np.linspace(-Ta, 125, 500)

	fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3.6))
	fig.subplots_adjust(top=0.92, bottom=0.155, left=0.075, right=0.985, hspace=0.2, wspace=0.2)

	for theta, ax in zip([0, 60], [ax1, ax2]):
		plot_failure_curve(ax, 0.0, theta, I1_list, gamma, beta, beta1, Ta)
		plot_failure_curve(ax, 0.0003, theta, I1_list, gamma, beta, beta1, Ta)
		plot_failure_curve(ax, 0.0006, theta, I1_list, gamma, beta, beta1, Ta)
		plot_failure_curve(ax, 0.0008, theta, I1_list, gamma, beta, beta1, Ta)
		plot_failure_curve(ax, 0.0016, theta, I1_list, gamma, beta, beta1, Ta)
		plot_dilatancy_boundary(ax, theta, I1_list, gamma, beta, beta1, Ta)
		ax.set_xlabel(r"$I_1$" + " (MPa)")
		ax.set_ylabel(r"$\sqrt{J_2}$" + " (MPa)")
		ax.grid(True)
		ax.set_ylim(0, 45)
		ax.set_title(r"$\theta =$" + str(theta) + r"$^o$")

	ax1.text(53, 2.6, r"$\alpha=0.0016$", fontsize=8)
	ax1.text(100.8, 8.4, r"$\alpha=0.0008$", fontsize=8)
	ax1.text(80.4, 15.6, r"$\alpha=0.0006$", fontsize=8)
	ax1.text(106, 24.4, r"$\alpha=0.0003$", fontsize=8)
	ax1.text(106, 30.6, r"$\alpha=0.0$", fontsize=8)
	ax1.text(-5, 16.2, "Dilatancy boundary", fontsize=9, fontname="serif")
	ax1.text(33.6, 26.3, "Failure boundary", fontsize=9, fontname="serif")
	ax1.plot([41.8, 21], [7.2, 15.1], "k-", linewidth=0.5)
	ax1.plot([73.5, 53.5], [19.8, 25.4], "k-", linewidth=0.5)

	ax2.text(52.6, 5.6, r"$\alpha=0.0016$", fontsize=8)
	ax2.text(102.9, 12, r"$\alpha=0.0008$", fontsize=8)
	ax2.text(76, 21.8, r"$\alpha=0.0006$", fontsize=8)
	ax2.text(100.7, 32.2, r"$\alpha=0.0003$", fontsize=8)
	ax2.text(100.7, 40.6, r"$\alpha=0.0$", fontsize=8)
	ax2.text(-7.6, 22.6, "Dilatancy boundary", fontsize=9, fontname="serif")
	ax2.text(10.1, 30.9, "Failure boundary", fontsize=9, fontname="serif")
	ax2.plot([39.2, 18.4], [12.0, 21.8], "k-", linewidth=0.5)
	ax2.plot([59.6, 30.1], [26.3, 30.4], "k-", linewidth=0.5)

	apply_white_theme(fig, [ax1, ax2], transparent=True)

	plt.show()


if __name__ == '__main__':
	main2()