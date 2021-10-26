using Plots; plotly

n=10
nt = 60
psi_partions_num = 20
x_res=n #number of cells in x dim
y_res=2
in_data=zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
length_domain = 0.18 #length of periodic element
height_domain = 0.04
dt=0.005
read!("PSP_C15_uniform_1_thin_c_0_10_new_abs_inf_vp_CLT",in_data)
t_space = LinRange(0,nt*dt,nt)
x_space = LinRange(0,length_domain,x_res)
y_space = LinRange(0,height_domain,y_res)

psi_spacing = LinRange(0,1.2,psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

phi_1_means = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt]
phi_2_means = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing ) for psi_1_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt]

x_integral_means = sum(phi_1_means, dims=2)[:,1,:] .* 1/x_res

x_integral_f_phi_1 = sum(in_data,dims=[2,4])[:,1,:,1,:].* 1/x_res

# surface(t_space,y_space,x_integral_means)
# plot(y_space,x_integral_means[:,nt])
# surface(x_space,y_space,phi_2_means[:,:,nt], ylabel="y")
# surface(t_space,psi_spacing,x_integral_f_phi_1[:,1,:])
plot(psi_spacing,x_integral_f_phi_1[:,1,nt], ylabel="Relative frequency", xlabel="ϕ",
title="pdf of ϕ at time 0.3 for N → ∞ with PSP mixing")
savefig("CLT_reaction_Ninf_PSP.png")
# plot(t_space,x_integral_means[1,:])