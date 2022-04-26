using Plots; GR.init()#plotlyjs()
import Statistics as st

data_folder=  "Data/paper/toy_prob_9/with_phi_gamma/F32/"#"Data/default/F32/"

image_path=["plots"]

image_folder=""
for folder in image_path
    global image_folder*=folder
    try
        run(`mkdir $image_folder`)
    catch e 
        isa(e,ProcessFailedException) || rethrow()
    end
    global image_folder*="/"
end

float_type=Float32

data_shape = zeros(Int,5)
# read!(data_folder*"array_shape",data_shape)
#  [psi_partions_num, psi_partions_num, y_res, x_res, nt+1]
data_shape = [50,50,5,15,114]
nt = data_shape[5]
psi_partions_num = data_shape[1]
phi_domain = [0,1.01]
Del_phi = (phi_domain[2]-phi_domain[1])/psi_partions_num
Del_phi_x_phi = (Del_phi)^2 
Del_phi_x_phi = (Del_phi)^2 
Del_phi_x_phi = (Del_phi)^2 
x_res=data_shape[4] #number of cells in x dim
y_res=data_shape[3]
in_data=zeros(float_type,psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
length_domain = 0.20 #length of periodic element
height_domain = 0.01
T=0.1 #total time
dt=T/(nt)
k = 0.938
C_0 = 2.1
omega_mean=5
B=(1.5*C_0)
bc_k=0.25
D=2*C_0*k/((B^2)*omega_mean)
do_flux = false
# read!(data_folder*"data",in_data)
read!(data_folder*"Inf_vp_CLT_w50.0",in_data)
t_space = LinRange(0,nt*dt,nt+1)
x_space = LinRange(0,length_domain,x_res)
y_space = LinRange(0,height_domain,y_res)

psi_spacing = LinRange(phi_domain[1],phi_domain[2],psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

phi_1_means = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]
phi_2_means = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing ) for psi_1_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]

x_integral_means_1 = sum(phi_1_means, dims=2)[:,1,:] .* 1/x_res
y_integral_means_1 = sum(phi_1_means, dims=1)[1,:,:] .* 1/y_res

x_integral_f_phi_1 = sum(in_data,dims=[2,4])[:,1,:,1,:].* 1/x_res
integral_f_phi_1 = sum(x_integral_f_phi_1,dims=2)[:,1,:].* 1/y_res
integral_mean_1 = sum(phi_1_means,dims=[1,2])[1,1,:].* 1/(y_res*x_res)

plt1=plot()
plt2=plot()
plt3=plot()

plot!(plt1,t_space,integral_mean_1,
title="Evolution of mean concentration over time",
ylabel="concentration",
xlabel="time",
label=""
)

contour!(plt2,t_space[:],psi_spacing,integral_f_phi_1[:,:],
levels=200, 
color= colormap("RdBu", logscale=true),
title="Evolution of whole domain PDF over time",
ylabel="concentration",
xlabel="time",
colorbartitle="Relative frequency",
fill=true,
line=false)

contour!(plt3,t_space[:],y_space,x_integral_means_1[:,:],
levels=200, 
# color= colormap("RdBu", logscale=true),
title="Evolution of mean concentration vs y over time",
ylabel="y",
xlabel="time",
colorbartitle="concentration",
fill=true,
line=false)

savefig(plt1,image_folder*"mean_over_time")
savefig(plt2,image_folder*"pdf_over_time")
savefig(plt3,image_folder*"mean_y_over_time")