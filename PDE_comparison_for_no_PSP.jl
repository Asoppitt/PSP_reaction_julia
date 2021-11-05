using OrdinaryDiffEq, ModelingToolkit, DiffEqOperators, Plots; plotlyjs()

k = 0.938
C_0 = 3 
omega_mean=100
B=(1.5*C_0)
bc_k=0.25
D=2*C_0*k/((B^2)*omega_mean)
height_domain =0.04

nknots = 200
h = height_domain/(nknots+1)
knots = range(0, step=h, length=nknots)
ord_deriv = 2
ord_approx = 2

const Δ = D*CenteredDifference(ord_deriv, ord_approx, h, nknots)
const bc = RobinBC((-bc_k,D,0.0),(0.0,D,0.0), h, ord_approx)

t0 = 0.0
t1 = 0.3
u0 = ones(size(knots)) #.*(knots.<0.5*height_domain)#u_analytic.(knots, t0)

step(u,p,t) = Δ*bc*u
prob = ODEProblem(step, u0, (t0, t1))
alg = KenCarp4()
sol = solve(prob, alg)

U = [sol.u[i][j] for i in 1:length(sol.t), j in 1:nknots]
surface(knots,sol.t,U,ylabel="t",
xlabel="y")