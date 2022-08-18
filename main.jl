using FractionalDiffEq, Plots

# Fractiona delay differential equations with constant delay
ϕ(x) = x == 0 ? 19.00001 : 19.0
f(t, y, ϕ) = 3.5*y*(1-ϕ/19)
h = 0.05; α = 0.97; τ = 0.8; T = 56
fddeprob = FDDEProblem(f, ϕ, α, τ, T)
V, y = solve(fddeprob, h, DelayPECE())
plot(y, V, xlabel="y(t)", ylabel="y(t-τ)")

# Fractional delay differential equations with multiple delays
using FractionalDiffEq, Plots
α = 0.95; ϕ(x) = 0.5
τ = [2, 2.6]
fun(t, y, ϕ1, ϕ2) = 2*ϕ1/(1+ϕ2^9.65)-y
prob = FDDEProblem(fun, ϕ, α, τ, 100)
delayed, y = solve(prob, 0.01, DelayPECE())

p1=plot(delayed[1, :], y)
p2=plot(delayed[2, :], y)
plot(p1, p2, layout=(1, 2))

# System of fractional differential equations
using FractionalDiffEq, Plots
α=[0.94, 0.94, 0.94]; ϕ=[0.2, 0, 0.5]; τ=0.009; T=1.4; h=0.001
function delaychen!(dy, y, ϕ, t)
	a=35; b=3; c=27
	dy[1] = a*(y[2]-ϕ[1])
	dy[2] = (c-a)*ϕ[1]-y[1]*y[3]+c*y[2]
	dy[3] = y[1]*y[2]-b*ϕ[3]
end
prob = FDDESystem(delaychen!, ϕ, α, τ, T)
sol=solve(prob, h, DelayABM())
plot(sol, title="Fractional Order Chen Delayed System")


# Matrix for fractional delay differential equations
using FractionalDiffEq, Plots

tspan = (0, 70); τ=3.1416; h=0.01; α=0.4
x0(t) = [sin(t)*cos(t); sin(t)*cos(t); cos(t)^2-sin(t)^2; cos(t)^2-sin(t)^2]
A=[0 0 1 0; 0 0 0 1; 0 -2 0 0; -2 0 0 0]
B=[0 0 0 0; 0 0 0 0; -2 0 0 0; 0 -2 0 0]
f=[0; 0; 0; 0]

prob = FDDEMatrixProblem(α, τ, A, B, f, x0, tspan)
sol=solve(prob, h, MatrixForm())
plot(sol[:, 1], sol[:, 3])