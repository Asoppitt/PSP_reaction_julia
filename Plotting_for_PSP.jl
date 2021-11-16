using Plots; plotlyjs()#GR.init()

n=10
nt = 60
psi_partions_num = 20
phi_domain = [0,1.2]
Del_phi = (phi_domain[2]-phi_domain[1])/psi_partions_num
Del_phi_x_phi = (Del_phi)^2 
x_res=n #number of cells in x dim
y_res=10*n
in_data=zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
in_data_flux=zeros(2,x_res, nt)
in_data_flux_vel=zeros(x_res, nt)
length_domain = 0.18 #length of periodic element
height_domain = 0.04
dt=0.005
k = 0.938
C_0 = 2.1
omega_mean=10
B=(1.5*C_0)
bc_k=0.25
D=2*C_0*k/((B^2)*omega_mean)
do_flux = true
filename = "Data/PSP_off_uniform_1_tallinreslonginspace_c_0_21_k_08_w_1_new_abs_80_vp_CLT"
read!(filename,in_data)
do_flux && read!(filename*"flux",in_data_flux)
do_flux && read!(filename*"flux_vel",in_data_flux_vel)
t_space = LinRange(0,nt*dt,nt)
x_space = LinRange(0,length_domain,x_res)
y_space = LinRange(0,height_domain,y_res)

psi_spacing = LinRange(0,1.2,psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

phi_1_means = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt]
phi_2_means = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing ) for psi_1_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt]

x_integral_means = sum(phi_1_means, dims=2)[:,1,:] .* 1/x_res

x_integral_f_phi_1 = sum(in_data,dims=[2,4])[:,1,:,1,:].* 1/x_res./Del_phi

flux_y0_vel = sum(in_data_flux_vel[:,:], dims=1)[1,:]
flux_y0_1_int = sum(in_data_flux[1,:,:], dims=1)[1,:].* 1/x_res
hopefully_constant  = (flux_y0_1_int)./(x_integral_means[1,:])


# do_flux && plot!(t_space,flux_y0_1_int)
do_flux && plot(t_space,flux_y0_vel)
# surface(t_space,y_space,x_integral_means,ylabel="y",
# xlabel="t",title="Particle method")
# plot(y_space,x_integral_means[:,nt])
# surface(x_space,y_space,phi_2_means[:,:,nt], ylabel="y")
# heatmap(t_space,psi_spacing,x_integral_f_phi_1[:,1,1:nt], levels=50,
# color= colormap("RdBu", logscale=true),
# title="Contour Plot of ϕ denisty over time for N → ∞",
# ylabel="ϕ",
# xlabel="time",
# colorbartitle="Relative frequency",)
# contour(t_space,psi_spacing,x_integral_f_phi_1[:,1,1:nt], levels=100, 
# color= colormap("RdBu", logscale=true),
# title="Contour Plot of ϕ denisty over time for N → ∞",
# ylabel="ϕ",
# xlabel="time",
# colorbartitle="Relative frequency",
# fill=true,
# line=true)
# plot(psi_spacing,x_integral_f_phi_1[:,1,20], ylabel="Relative frequency", xlabel="ϕ",
# title="pdf of ϕ at time 0.3 for N = 1 with PSP mixing",label="2")
# savefig("CLT_reaction_N1_PSP.png")
# plot(t_space,x_integral_means[1,:])