using Plots
using Roots

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
function radial_wfc_numerov_at0(E, r0, n, l, Z)
	ur = zeros(size(r0))
	fn = zeros(size(r0))
	ur[end] = 0.0
	ur[end-1] = 0.001

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
	# now extrapolate the wavefunction to u(0)
	# the first derivative equals at ur[0]:(u0 - ur[0]) / (0 - r0[0]) = (ur[1] - ur[0]) / (r0[1] - r0[0])
	u0 = ur[1] + (ur[2] - ur[1]) * (0 - r0[1]) / dr
	return u0
end

function main()
	Z = 1.0
	n = 3.0
	l = 2
	r0 = range(1e-6, stop = 30, length = 3000) 
	E_lower = E_upper = -0.14 # In practice, the guess energy should be smaller than the smallest potential energy -Z^2/n^2.
	dE = 0.01
	println("The guess initial energy is:E_lower=$E_lower,E_upper=$E_upper")
	u1 = radial_wfc_numerov_at0(E_lower, r0, n, l, Z)
	println("shooting:")
	while true
		E_upper += dE
		println("E_lower=$E_lower,E_upper=$E_upper")
		u2 = radial_wfc_numerov_at0(E_upper, r0, n, l, Z)
		if u1 * u2 < 0
			break
		end
	end
	println("The energy is locating on [$E_lower,$E_upper]")
	E = find_zero((E) -> radial_wfc_numerov_at0(E, r0, n, l, Z), (E_lower, E_upper), Roots.Brent())
	println("The energy is:E = $E")
end

main()




