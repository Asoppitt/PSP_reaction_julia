using OrdinaryDiffEq, ModelingToolkit, DiffEqOperators, Plots; plotlyjs()

k = 0.938
C_0 = 2.1
omega_mean=10.0
B=(1.5*C_0)
bc_k=0.3
D=2*C_0*k/((B^2)*omega_mean)
height_domain =0.01
length_domain =0.20
dt=0.001
nt=200
T=dt*nt

nknots = 200
h = height_domain/(nknots+1)
knots = range(0, step=h, length=nknots)
ord_deriv = 2
ord_approx = 2

const Δ = (D)*CenteredDifference(ord_deriv, ord_approx, h, nknots)
const bc = RobinBC((-bc_k,D,0.0),(0.0,D,0.0), h, ord_approx)

t0 = 0.0
t1 = T
u0 = ones(size(knots)) #.*(knots.<0.5*height_domain)#u_analytic.(knots, t0)

step(u,p,t) = Δ*bc*u
prob = ODEProblem(step, u0, (t0, t1),reltol=1e-8, abstol=1e-8)
alg = KenCarp4()
sol = solve(prob, alg)

U = [sol.u[i][j] for j in 1:nknots, i in 1:length(sol.t)]
# plt1 = surface(sol.t,knots,U,ylabel="y",
# xlabel="t",title="MOL solver")
# display(plt1)

plt2= plot(sol.t,U[1,:])
# plot!(sol.t,U[end,:])
# plot!(sol.t,sum(U[:,:],dims=1)[1,:]*1/nknots)
display(plt2)

# plt3= plot((1:nknots).*h,U[:,end])
# display(plt3)