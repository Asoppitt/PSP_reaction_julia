using OrdinaryDiffEq, ModelingToolkit, DiffEqOperators, Plots

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
t0 = 0.0
t1 = T
ord_deriv = 2
ord_approx = 2
big_run = 2^10
h_big = height_domain/(big_run-1)
reltol=1e-12
abstol=1e-12
alg = RadauIIA5()

# for nknots=[big_run]
#     h=h_big
#     knots = range(0, step=h, length=nknots)

#     Δ = (D)*CenteredDifference(ord_deriv, ord_approx, h, nknots)
#     bc = RobinBC((-bc_k,D,0.0),(0.0,D,0.0), h, ord_approx)

#     u0 = ones(size(knots)) #.*(knots.<0.5*height_domain)#u_analytic.(knots, t0)

#     step(u,p,t) = Δ*bc*u
#     prob = ODEProblem(step, u0, (t0, t1),reltol=1e-12, abstol=1e-12)
#     global sol_large_nknots = solve(prob, alg)
# end
# integration_weights_1d = sol_large_nknots.t[2:end]-sol_large_nknots.t[1:end-1]
# integration_weights=[h_big*integration_weights_1d[k] for j=1:big_run, k=1:length(integration_weights_1d)]

nknots_array=2 .^(2:7)
error_num=zeros(length(nknots_array))
# for (i,nknots) in enumerate(nknots_array)
#     h = height_domain/(nknots-1)
#     knots = range(0, step=h, length=nknots)

#     Δ = (D)*CenteredDifference(ord_deriv, ord_approx, h, nknots)
#     bc = RobinBC((-bc_k,D,0.0),(0.0,D,0.0), h, ord_approx)

#     u0 = ones(size(knots)) #.*(knots.<0.5*height_domain)#u_analytic.(knots, t0)

#     step(u,p,t) = Δ*bc*u
#     prob = ODEProblem(step, u0, (t0, t1),reltol=reltol, abstol=abstol)
#     sol = solve(prob, alg)
#     sol_at_points=sol.(sol_large_nknots.t[2:end])
#     sol_at_points=[sol_at_points[k][Int(ceil(j*nknots//big_run))] for j=1:big_run, k=1:length(sol_large_nknots.t)-1]

#     error_num[i] = sum(integration_weights.*(sol_at_points-sol_large_nknots[:,2:end]).^2)
# end

read!("conv_study_error", error_num)

plt=plot(nknots_array,error_num)
display(plt)
plt2=plot(nknots_array,log10.(error_num))
display(plt2)