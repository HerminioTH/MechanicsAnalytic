import sympy as sy
import numpy as np
import matplotlib.pyplot as plt

MPa = 1e6

# def compute_J2_for_F_equal_0(I1, theta_deg, gamma=0.11, alpha=0.0, n=3, m=-0.5, beta=0.995, beta1=4.8e-9, Ta=5.4*MPa):
# 	theta_rad = np.radians(theta_deg)
# 	I1_star = I1 + Ta
# 	J2 = (-alpha*I1_star**n + gamma*I1_star**2)*(np.exp(beta1*I1_star) - beta*np.cos(3*theta_rad))**m
# 	print(J2)
# 	np.clip(J2, a_min=0.0, a_max=None, out=J2)
# 	return np.sqrt(J2)
# 	# return J2

def compute_J2_for_F_equal_0(I1, theta_deg, alpha=0.0, gamma=0.11, n=3, m=-0.5, beta=0.995, beta1=4.8e-9, Ta=5.4):
	theta_rad = np.radians(theta_deg)
	I1_star = I1 + Ta
	J2 = (-alpha*I1_star**n + gamma*I1_star**2)*(np.exp(beta1*I1_star) - beta*np.cos(3*theta_rad))**m
	np.clip(J2, a_min=0.0, a_max=None, out=J2)
	return np.sqrt(J2)

def dilatancy_boundary(I1, theta_deg, gamma=0.11, n=3, m=-0.5, beta=0.995, beta1=4.8e-9, Ta=5.4):
	theta_rad = np.radians(theta_deg)
	I1_star = I1 + Ta
	J2 = (1 - 2/n)*gamma*(I1_star**2)*(np.exp(beta1*I1_star) - beta*np.cos(3*theta_rad))**m
	return np.sqrt(J2)

def plot_failure_curve(ax, alpha, theta_deg, I1_list, gamma=0.11, beta=-0.995, beta1=4.8e-3, Ta=5.4/3):
	J2_sqr = compute_J2_for_F_equal_0(I1_list, alpha=alpha, theta_deg=theta_deg, gamma=gamma, beta=beta, beta1=beta1, Ta=Ta)
	ax.plot(I1_list, J2_sqr, "-", color="steelblue")

def plot_dilatancy_boundary(ax, theta_deg, I1_list, gamma=0.11, beta=-0.995, beta1=4.8e-3, Ta=5.4/3):
	J2_dil = dilatancy_boundary(I1_list, theta_deg=theta_deg, gamma=gamma, beta=beta, beta1=beta1, Ta=Ta)
	ax.plot(I1_list, J2_dil, "--", color="lightcoral")

def main():
	# beta = -0.995
	# beta1 = 4.8e-3
	# gamma = 0.11
	# Ta = 5.4

	beta = -0.995
	beta1 = 0.00479
	gamma = 0.0945
	Ta = 5.4#1.79

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


if __name__ == '__main__':
	main()