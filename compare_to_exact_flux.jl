using Plots; plotlyjs()#GR.init()
import Statistics as st

n=5
nt = 60
psi_partions_num = 100
phi_domain = [0,1.2]
Del_phi = (phi_domain[2]-phi_domain[1])/psi_partions_num
Del_phi_x_phi = (Del_phi)^2 
x_res=20*n #number of cells in x dim
y_res=n
in_data=zeros(psi_partions_num, psi_partions_num, 1, x_res, nt)
in_data_f=zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
length_domain = 0.20 #length of periodic element
height_domain = 0.01
dt=0.002
k = 0.938
C_0 = 2.1
omega_mean=10.0
u_mean=1.0
B=(1.5*C_0)
bc_k=0.25
D=2*C_0*k/((B^2)*omega_mean)
base_filename="Data/PSP_on_uniform_1_vlong_c_0_21_k_09_w_10_new_abs_Inf_vp_psi100_record_flux_u006"

ts=[i*dt for i=0:nt-1]
xs = (1:x_res).*(length_domain/x_res)

psi_spacing = LinRange(phi_domain[1],phi_domain[2],psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

N_trials=2

flux_means = zeros(x_res,nt,N_trials)
phi_means = zeros(x_res,nt,N_trials)

for i =1:N_trials
    filename = base_filename*"/run"*string(i)
    read!(filename,in_data)
    read!(filename*"f_phi",in_data_f)
    flux_dist = in_data[:,:,1,:,:] 
    phi_dist = in_data_f[:,:,1,:,:] 
    flux_means[:,:,i] = [sum([sum(flux_dist[:,psi_2_index,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for j in 1:x_res, t in 1:nt]
    phi_means[:,:,i] = [sum([sum(phi_dist[:,psi_2_index,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for j in 1:x_res, t in 1:nt]
end

for t=[60]#10:10:60
    plt = plot()
    plt2 = plot()
    plt3 = plot()
    for i=1:N_trials
        plot!(plt, xs, flux_means[:,t,i])
        plot!(plt2, xs, phi_means[:,t,i])
        plot!(plt3, xs, phi_means[:,t,i].-flux_means[:,t,i])
    end
    mean_flux = st.mean(flux_means[:,t,:],dims=2)[:,1]
    a_f=50
    c_f=0.077
    b_f = (mean_flux[1]-c_f)*exp( a_f*xs[1] )
    exp_approx_flux(x) = exp(-a_f*x)*b_f+c_f
    plot!(plt, xs, mean_flux,ls=:dash,lw=3,lc=:red)
    plot!(plt, xs, exp_approx_flux.(xs),ls=:dash,lw=3,lc=:blue)
    title!(plt, "Just Flux")
    mean_phi = st.mean(phi_means[:,t,:],dims=2)[:,1]
    a_p=50
    c_p=0.11
    b_p=(mean_phi[1]-c_p)*exp( a_p*xs[1] )
    exp_approx_phi(x) = exp(-a_p*x)*b_p+c_p
    plot!(plt2, xs, mean_phi,ls=:dash,lw=3,lc=:red)
    plot!(plt2, xs, exp_approx_phi.(xs),ls=:dash,lw=3,lc=:blue)
    title!(plt2, "Just Phi")
    mean_diff = st.mean(phi_means[:,t,:].-flux_means[:,t,:],dims=2)[:,1]
    a_d=50
    c_d=0.034
    b_d=(mean_diff[1]-c_d)*exp( a_d*xs[1] )
    exp_approx_diff(x) = exp(-a_d*x)*b_d+c_d
    plot!(plt3, xs, mean_diff,ls=:dash,lw=3,lc=:red)
    plot!(plt3, xs, exp_approx_diff.(xs),ls=:dash,lw=3,lc=:blue)
    title!(plt3, "Difference")
    display(plt)
    display(plt2)
    display(plt3)
    println(b_f/b_p, ' ', bc_k/D)
    println(b_d/b_p, ' ', bc_k/D)
    println(b_d/b_f, ' ', bc_k/D)
end

while true
    println("waiting here")
    sleep(5)
end
# filename = base_filename*"/run"*string(1)
# read!(filename,in_data)
# f_phi_1 = sum(in_data[:,:,1,:,:], dims=2)[:,1,:,:]
# plt3 = surface((1:x_res).*(length_domain/x_res),psi_spacing,[f_phi_1[i,j,1]>0 ? log(f_phi_1[i,j,1]) : -10.0 for i=1:psi_partions_num, j=1:x_res], levels=100, 
# color= colormap("RdBu", logscale=true),
# fill=true,
# line=false)
# display(plt3)

# plt4 = surface(ts,psi_spacing,[st.mean(f_phi_1[i,:,j])>0 ? log(st.mean(f_phi_1[i,:,j])) : -11.0 for i=1:psi_partions_num, j=1:nt], levels=100, 
# color= colormap("RdBu", logscale=true),
# fill=true,
# line=false)
# display(plt4)

# filename = base_filename*"/run"*string(1)*"f_phi"
# read!(filename,in_data_f)
# f_phi_1 = sum(in_data_f[:,:,1,:,:], dims=2)[:,1,:,:]
# plt5 = surface((1:x_res).*(length_domain/x_res),psi_spacing,[f_phi_1[i,j,40]>0 ? log(f_phi_1[i,j,40]) : -10.0 for i=1:psi_partions_num, j=1:x_res], levels=100, 
# color= colormap("RdBu", logscale=true),
# fill=true,
# line=false)
# display(plt5)

# plt7 = surface(ts,psi_spacing,[f_phi_1[i,1,j]>0 ? log(f_phi_1[i,1,j]) : -11.0 for i=1:psi_partions_num, j=1:nt],
# levels=100, 
# color= colormap("RdBu", logscale=true),
# fill=true,
# line=false)
# display(plt7)