"""
'module HydrogenIonContinuum'
... a model of that contains methods for computing one-electron continuum,matrix elements,etc..
"""

using SpecialFunctions
using HypergeometricFunctions
using Plots
struct QuantumNumber
      energy::Float64
      angularL::Int64
end
function radilaOrbit(qm::QuantumNumber, Z::Int64, r::Float64)
      E = qm.energy
      l = qm.angularL
      k = sqrt(2 * E)
      value = 2 * k * exp(pi * Z / (2 * k)) / factorial(2 * l + 1) * (2 * k * r)^(l)*exp(-k*r*im) * abs(gamma(l + 1 + Z / k * im))
      value = value * HypergeometricFunctions.M(l + 1 + Z / k * im, 2 * l + 2, 2 * k * r * im)
      return (value)
end

#s,p,d states
E=1
l=1
qm = QuantumNumber(E, l)
r_size = 0:0.01:20
Z = 1
R = []
for r in r_size
    push!(R, radilaOrbit(qm, Z, r))
end

R_r=real(R)
plot!(r_size,R_r)
