function hydrogen_potential(r, Z, n, l)
	# Hydrogen atom potential with angular momentum term
	return -Z^2 / n^2 - (-2.0 / r) - l * (l + 1) / (r^2)
end
function simpson_rule_squared(x, f)
	n = length(x)
	h = (x[end] - x[1]) / (n - 1)  # 计算小区间的宽度

	integral_sum = f[1]^2 + f[end]^2  # 加上首尾两个点的平方值

	# 奇数项的加权和
	for i in 1:2:n-2
		integral_sum += 4 * f[i+1]^2
	end

	# 偶数项的加权和
	for i in 2:2:n-3
		integral_sum += 2 * f[i+1]^2
	end

	integral_sum *= h / 3.0  # 乘以小区间宽度的权重
	return integral_sum
end

function numerov_solver(V, r_values, h)
	# Initialize arrays to store wavefunction values
	N = length(r_values)
	psi = zeros(N)

	# Initial conditions for Numerov method
	psi[1] = 0.0
	psi[2] = 0.001

	# Numerov algorithm - forward
	for n in 3:N
		h_sq = h^2
		g_n_minus_1 = V(r_values[n-2], l)
		g_n = V(r_values[n-1], l)
		g_n_plus_1 = V(r_values[n], l)

		A = 1.0 + (h_sq / 12) * g_n_plus_1
		B = 2.0 * (1.0 - (5 * h_sq / 12) * g_n)
		C = 1.0 + (h_sq / 12) * g_n_minus_1

		psi[n] = (B * psi[n-1] - C * psi[n-2]) / A
	end
    # normalization using Simpson's rule
	normalization = sqrt(simpson_rule_squared(r_values, psi))
	psi ./= normalization
	return psi
end

# Example parameters
r_min = 1E-11
r_max = 40.0
num_points = 500000  #grid set: psi_10 is wrong
h_step = (r_max - r_min) / (num_points - 1)

r_values = range(r_min, stop = r_max, length = num_points)

# Angular momentum quantum number (adjust as needed)
Z = 1
n = 3
l = 2

# Solve Schrödinger equation for hydrogen atom using Numerov method
wavefunction = numerov_solver((r, l) -> hydrogen_potential(r, Z, n, l), r_values, h_step)

# Plot the result
using Plots
plot(r_values, wavefunction, label = "Wavefunction", xlabel = "r", ylabel = "ψ_10(r)", legend = :bottomright)
