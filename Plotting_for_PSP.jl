using Plots; GR.init()#plotlyjs()
using LinearAlgebra
using LinearRegression
import Statistics as st

data_folder= "Data/test_DNS_Params/F32/"

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
read!(data_folder*"array_shape",data_shape)
#  [psi_partions_num, psi_partions_num, y_res, x_res, nt+1]
# data_shape = [50,50,5,15,114]
nt = data_shape[5]-1
psi_partions_num = data_shape[1]
phi_domain = [0,1.01]
x_res=data_shape[4] #number of cells in x dim
y_res=data_shape[3]
in_data=zeros(float_type,psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
length_domain = 0.2 #length of periodic element
height_domain = 0.01
T=0.05 #total time
dt=T/(nt)
read!(data_folder*"data",in_data)
t_space = LinRange(0,nt*dt,nt+1)
x_space_edges = LinRange(0,length_domain,x_res+1)
x_space = 0.5*(x_space_edges[2:end]+x_space_edges[1:end-1])
y_space_edges = LinRange(0,height_domain,y_res+1)
y_space = 0.5*(y_space_edges[2:end]+y_space_edges[1:end-1])

psi_spacing = LinRange(phi_domain[1],phi_domain[2],psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

phi_1_means = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]
# phi_2_means = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing ) for psi_1_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]

phi_1_2ndmom = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing.^2 ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]
# phi_2_2ndmom = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing.^2 ) for psi_1_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]

x_integral_means_1 = sum(phi_1_means, dims=2)[:,1,:] .* 1/x_res
# y_integral_means_1 = sum(phi_1_means, dims=1)[1,:,:] .* 1/y_res

phi_1_means_grad = [ind ? (phi_1_means[i+1,j,t]-phi_1_means[i-1,j,t])*2*(height_domain/y_res) :
                    (phi_1_means[i,j+1,t]-phi_1_means[i,j-1,t])*2*(length_domain/x_res) 
                    for ind in [true,false], i in 2:y_res-1, j in 2:x_res-1, t in 1:nt+1]

chi_grad=sum([norm(phi_1_means_grad[:,i,j,t]) for i in 1:y_res-2, j in 1:x_res-2, t in 1:nt+1],dims=[1,2])[1,1,:]

x_integral_f_phi_1 = sum(in_data,dims=[2,4])[:,1,:,1,:].* 1/x_res
integral_f_phi_1 = sum(x_integral_f_phi_1,dims=2)[:,1,:].* 1/y_res
integral_mean_1 = sum(phi_1_means,dims=[1,2])[1,1,:].* 1/(y_res*x_res)
integral_2ndmom_1 = sum(phi_1_2ndmom,dims=[1,2])[1,1,:].* 1/(y_res*x_res)*length_domain*height_domain

plt1=plot()
plt1_1=plot()
plt1_12=plot()
plt1_2=plot()
plt2=plot()
plt3=plot()
plt4=plot()
plt5=plot()


D=2

chi_exact_t = integral_mean_1[1]^2*height_domain/(8*sqrt(2*pi*D))*t_space.^(-1.5)

plot!(plt1,t_space,chi_grad,
title="chi grad",
ylabel="concentration",
xlabel="time",
label=""
)

# dmdt=-0.5*(integral_mean_1[3:end]-integral_mean_1[1:end-2])/dt

max_grad=argmax(chi_grad)

end_time_grad=nt
try
    end_time_grad=minimum(findall(chi_grad[max_grad:end].<0))-1
catch e
    isa(e,ArgumentError) || rethrow()
end

power_g, constant_g = linregress(log.(t_space[max_grad:end_time_grad]),log.(chi_grad[max_grad:end_time_grad])).coeffs

println("power_grad=",power_g)

scatter!(plt1_1,log.(t_space[max_grad:end_time_grad]),log.(chi_grad[max_grad:end_time_grad]),
title="log log chi_grad",
ylabel="reaction rate",
xlabel="time",
label=""
)

plot!(plt1_1,log.(t_space[max_grad:end_time_grad]),power_g.*log.(t_space[max_grad:end_time_grad]).+constant_g,
label=""
)

plot!(plt1_1,log.(t_space[max_grad+1:end_time_grad]),log.(chi_exact_t[max_grad+1:end_time_grad]))

# plot!(plt1_12,t_space[2:end-1],dmdt,
# title="Comparison of scalar disipation rate and reaction rate",
# ylabel="reaction rate",
# xlabel="time",
# label=""
# )


chi_t=-0.5*(integral_2ndmom_1[3:end]-integral_2ndmom_1[1:end-2])/dt

max_m=argmax(chi_t)

plot!(plt1_12,t_space[2:end-1],chi_t,label="χ")
# plot!(plt1_12,t_space[max_m:end],chi_exact_t[max_m:end],label="χ theoretical")

end_time=nt-1
try
    end_time=minimum(findall(chi_t[max_m:end].<0))-1
catch e
    isa(e,ArgumentError) || rethrow()
end

power, constant = linregress(log.(t_space[max_m:end_time]),log.(chi_t[max_m:end_time])).coeffs

println("power=",power)

scatter!(plt1_2,log.(t_space[max_m+1:end_time]),log.(chi_t[max_m:end_time-1]),
title="Evolution of scalar disipation rate over time",
ylabel="χ",
xlabel="time",
label=""
)

plot!(plt1_2,log.(t_space[max_m+1:end_time]),log.(t_space[max_m+1:end_time]).*power.+constant,label="")

plot!(plt1_2,log.(t_space[max_m+1:end_time]),log.(chi_exact_t[max_m+1:end_time]))

heatmap!(plt2,t_space[:],psi_spacing,integral_f_phi_1[:,:],
levels=200, 
color= colormap("RdBu", logscale=true),
title="Evolution of whole domain PDF over time",
ylabel="concentration",
xlabel="time",
colorbartitle="Relative frequency",
fill=true,
line=false)

# heatmap!(plt3,t_space[:],x_space,y_integral_means_1[:,:],
# levels=200, 
# # color= colormap("RdBu", logscale=true),
# title="Evolution of mean concentration vs y over time",
# ylabel="y",
# xlabel="time",
# colorbartitle="concentration",
# fill=true,
# line=false)

heatmap!(plt4,x_space[:],y_space,phi_1_means[:,:,1],#change the final index in phi_1_means to choose time plotted
levels=200, 
# color= colormap("RdBu", logscale=true),
title="Mean ϕ1 vs space at the intial time",
ylabel="y",
xlabel="x",
colorbartitle="concentration",
fill=true,
line=false)

heatmap!(plt5,x_space[:],y_space,phi_1_means[:,:,end],#change the final index in phi_1_means to choose time plotted
levels=200, 
# color= colormap("RdBu", logscale=true),
title="Mean ϕ1 vs space at the final time",
ylabel="y",
xlabel="x",
colorbartitle="concentration",
fill=true,
line=false)

for y in y_space_edges
    plot!(plt4, [0,length_domain],[y,y],label="",color="cyan")
end
for x in x_space_edges
    plot!(plt4, [x,x], [0,height_domain],label="",color="cyan")
end

# savefig(plt1,image_folder*"mean_over_time")
# savefig(plt1_1,image_folder*"rate_over_time")
# savefig(plt1_12,image_folder*"rvchi_over_time")
# savefig(plt1_2,image_folder*"chi_over_time")
# savefig(plt2,image_folder*"pdf_over_time")
# savefig(plt3,image_folder*"mean_y_over_time")
# savefig(plt4,image_folder*"initial_time_concentarions")
# savefig(plt5,image_folder*"final_time_concentarions")

# uncommentfor display. If Plots backend is set to plotlyjs these are interactive
display(plt1)
display(plt1_1)
display(plt1_2)
display(plt1_12)
display(plt2)
display(plt3)
display(plt4)
display(plt5)
# #used for holding host process for plots open - not needed for vscode operation
# print(stdout, "press enter to close plots")
# read(stdin, 1)