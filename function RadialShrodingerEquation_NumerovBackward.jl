using Plots


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

#=
y''(x) = - g(x)y(x)
Numerov algorithm

				  [12 - 10f(n)]*y(n) - y(n-1)*f(n-1)
		y(n+1) = ------------------------------------
							   f(n+1)

	where

		f(n) = 1 + (h**2 / 12)*g(n)

		g(n) = [E + (2*Z / x) - l*(l+1) / x**2]
=#
function radial_wfc_numerov(r0, n, l, Z, du = 0.001)
	ur = zeros(size(r0))
	fn = zeros(size(r0))

	E = -Z^2 / n^2
	ur[end] = 0.0
	ur[end-1] = du

	dr = r0[2] - r0[1]
	h12 = dr^2 / 12.0

	gn = (E .+ 2 * Z ./ r0 - l * (l + 1) ./ r0 .^ 2)
	fn = 1.0 .+ h12 * gn
	#backward
	for ii in reverse(1:length(r0)-2)
		ur[ii] = (12 - 10 * fn[ii+1]) * ur[ii+1] -
				 ur[ii+2] * fn[ii+2]
		ur[ii] /= fn[ii]
	end

	# normalization using Simpson's rule
	normalization = sqrt(simpson_rule_squared(r0, ur))
	ur ./= normalization

	return ur
end

# Define the radial grid
r0 = range(1e-10, stop = 40, length = 1000) #Note that Numerov’s method only support linear mesh, as can be seen in the method section.
# Set other parameters (use Float64 for Z, n, l)
Z = 1.0
n = 3.0
l = 0
# Backward integration from u(\infty) with Numerov method
ur3 = radial_wfc_numerov(collect(r0), n, l, Z)
plot(r0, ur3, label = "Wavefunction", xlabel = "r", ylabel = "ψ(r)", legend = :bottomright)
