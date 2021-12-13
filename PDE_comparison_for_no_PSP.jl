using OrdinaryDiffEq, ModelingToolkit, DiffEqOperators, Plots; plotlyjs()

k = 0.938
C_0 = 2.1
omega_mean=1.0
B=(1.5*C_0)
bc_k=1.0
D=0.1#2*C_0*k/((B^2)*omega_mean)
height_domain =0.04
length_domain =0.10
v_mean = 1.0

nknots = 200
h = height_domain/(nknots+1)
knots = range(0, step=h, length=nknots)
ord_deriv = 2
ord_approx = 2

const Δ = (D/v_mean)*CenteredDifference(ord_deriv, ord_approx, h, nknots)
const bc = RobinBC((-bc_k,D,0.0),(0.0,D,0.0), h, ord_approx)

x0 = 0.0
x1 = length_domain
u0 = ones(size(knots)) #.*(knots.<0.5*height_domain)#u_analytic.(knots, t0)

step(u,p,t) = Δ*bc*u
prob = ODEProblem(step, u0, (x0, x1),reltol=1e-8, abstol=1e-8)
alg = KenCarp4()
sol = solve(prob, alg)

U = [sol.u[i][j] for j in 1:nknots, i in 1:length(sol.t)]
plt1 = surface(sol.t,knots,U,ylabel="y",
xlabel="x",title="MOL solver")
display(plt1)

cubic_approx(x) = 1/x^(1/3)*(sol.t[10]^(1/3)*bc_k/D.*U[1,10])
a=25
exp_approx(x) = exp(-a*x)*(exp(sol.t[1])*bc_k/D.*U[1,1])
plt2= plot(sol.t,bc_k/D.*U[1,:])
plot!(plt2, sol.t[10:end],cubic_approx.(sol.t[10:end]))
plot!(plt2, sol.t,exp_approx.(sol.t))
display(plt2)
