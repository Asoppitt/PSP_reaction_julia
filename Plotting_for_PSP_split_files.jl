using Plots; GR.init()#plotlyjs()
using LinearAlgebra, Printf
using LinearRegression
import Statistics as st

data_folder= "Data/test_2_2/1500_part/F32/"

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

total_shape_info = [0,0,0,false,false]

try
    read!(data_folder*"total_shape",total_shape_info)
catch
    read!(data_folder*"total_shape",total_shape_info[1:4])
end

dim_phi=0.1;

record_moments=Bool(total_shape_info[5])
# record_moments=false

#get the number of chunks, chunk length and if final chunck is full length
n_chunks=total_shape_info[1]-total_shape_info[4]
nt=total_shape_info[3]
#get psi_partions_num and x/y _res from the first chunck
filename_chunck1=data_folder*"1_"*string(total_shape_info[2])
data_shape = zeros(Int,5)
read!(filename_chunck1*"array_shape",data_shape)

record_moments||(integral_f_phi_1=zeros(float_type,data_shape[[1]]...,nt))
integral_dist_1 = zeros(float_type, data_shape[1],nt)
mean_field_1=zeros(float_type,data_shape[[3,4]]...,nt)
mean_field_2=zeros(float_type,data_shape[[3,4]]...,nt)
prod_field=zeros(float_type,data_shape[[3,4]]...,nt)
record_moments&&(mean_1=zeros(float_type,nt);mean_2=zeros(float_type,nt))
record_moments&&(integral_2ndmom_1=zeros(float_type,nt);integral_2ndmom_2=zeros(float_type,nt))

psi_partions_num_1=data_shape[1]
psi_partions_num_2=data_shape[2]
x_res=data_shape[4] #number of cells in x dim
y_res=data_shape[3]
length_domain = 5e-3 #length of periodic element
height_domain = length_domain./240
T=250000*7.63364e-09 #total time
phi_domain=[0,10.01]
dt=T/100#nt#T/(nt)
t_space = LinRange(0,nt*dt,nt+1)
x_space_edges = LinRange(0,length_domain,x_res+1)
x_space = 0.5*(x_space_edges[2:end]+x_space_edges[1:end-1])
y_space_edges = LinRange(0,height_domain,y_res+1)
y_space = 0.5*(y_space_edges[2:end]+y_space_edges[1:end-1])

psi_spacing_1 = LinRange(phi_domain[1],phi_domain[2],psi_partions_num_1+1) #defines boundaries for cells in psi space
psi_spacing_1 = 0.5 .* (psi_spacing_1[1:psi_partions_num_1]+psi_spacing_1[2:psi_partions_num_1+1])#defining centres of bins

#use when processing binary data, gives more accurate means
psi_spacing_1 = Array(psi_spacing_1)

psi_spacing_2 = LinRange(phi_domain[1],phi_domain[2],psi_partions_num_2+1) #defines boundaries for cells in psi space
psi_spacing_2 = 0.5 .* (psi_spacing_2[1:psi_partions_num_2]+psi_spacing_2[2:psi_partions_num_2+1])#defining centres of bins

#use when processing binary data, gives more accurate means
psi_spacing_2 = Array(psi_spacing_2)
# psi_spacing[1] = 0
# psi_spacing[end] = 1

# target_cell=(1,110)

# single_cell_dist_t=zeros(float_type, psi_partions_num, nt)

#careful, this could be larger than ram
for chunk in 0:n_chunks-1
    filename=data_folder*string(chunk*total_shape_info[2]+1)*'_'*string((chunk+1)*total_shape_info[2])
    read!(filename*"array_shape",data_shape)
    in_data=zeros(float_type,data_shape...)
    read!(filename*"data",in_data)
    integral_dist_1_chunk = sum(in_data,dims=[2,3,4])[:,1,1,1,:].*1/(x_res*y_res)
    mean_field_1_chunk = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing_1 ) for psi_2_index in 1:psi_partions_num_2]) for i in 1:y_res, j in 1:x_res, t in 1:data_shape[5]]
    # mean_field_2_chunk = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing_2 ) for psi_1_index in 1:psi_partions_num_1]) for i in 1:y_res, j in 1:x_res, t in 1:data_shape[5]]
    mean_field_1[:,:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = mean_field_1_chunk
    # mean_field_2[:,:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = mean_field_2_chunk
    integral_dist_1[:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = integral_dist_1_chunk
    # prod_field_chunk = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing_1.*psi_spacing_2[psi_2_index] ) for psi_2_index in 1:psi_partions_num_2]) for i in 1:y_res, j in 1:x_res, t in 1:data_shape[5]]
    # prod_field[:,:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = prod_field_chunk
    if record_moments
        means_chunk = zeros(float_type,2,data_shape[5])
        mom2s_chunk = zeros(float_type,2,data_shape[5])
        read!(filename*"mean",means_chunk)
        read!(filename*"2nd_moment",mom2s_chunk)
        mean_1[(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = means_chunk[1,:]./dim_phi
        mean_2[(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = means_chunk[2,:]./dim_phi
        integral_2ndmom_1[(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = mom2s_chunk[1,:]./dim_phi^2
        integral_2ndmom_2[(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = mom2s_chunk[2,:]./dim_phi^2
    else
        integral_f_phi_1_chunk = sum(in_data,dims=[2,3,4])[:,1,1,1,:].* 1/y_res*1/x_res
        integral_f_phi_1[:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = integral_f_phi_1_chunk
    end
end

if total_shape_info[4]==1
    local chunk=n_chunks
    local filename=data_folder*string((n_chunks)*total_shape_info[2]+1)*'_'*string(nt)
    read!(filename*"array_shape",data_shape)   
    local in_data=zeros(float_type,data_shape...)
    read!(filename*"data",in_data)
    local integral_dist_1_chunk = sum(in_data,dims=[2,3,4])[:,1,1,1,:].*1/(x_res*y_res)
    local mean_field_1_chunk = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing_1 ) for psi_2_index in 1:psi_partions_num_2]) for i in 1:y_res, j in 1:x_res, t in 1:data_shape[5]]
    # local mean_field_2_chunk = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing_2 ) for psi_1_index in 1:psi_partions_num_1]) for i in 1:y_res, j in 1:x_res, t in 1:data_shape[5]]
    mean_field_1[:,:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = mean_field_1_chunk
    # mean_field_2[:,:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = mean_field_2_chunk
    integral_dist_1[:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = integral_dist_1_chunk

    # local prod_field_chunk = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing_1.*psi_spacing_2[psi_2_index] ) for psi_2_index in 1:psi_partions_num_2]) for i in 1:y_res, j in 1:x_res, t in 1:data_shape[5]]
    # prod_field[:,:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = prod_field_chunk
    if record_moments
        local means_chunk = zeros(float_type,2,data_shape[5])
        local mom2s_chunk = zeros(float_type,2,data_shape[5])
        read!(filename*"mean",means_chunk)
        read!(filename*"2nd_moment",mom2s_chunk)
        mean_1[(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = means_chunk[1,:]./dim_phi
        mean_2[(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = means_chunk[2,:]./dim_phi
        integral_2ndmom_1[(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = mom2s_chunk[1,:]./dim_phi^2
        integral_2ndmom_2[(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = mom2s_chunk[2,:]./dim_phi^2
    else
        local integral_f_phi_1_chunk = sum(in_data,dims=[2,3,4])[:,1,1,1,:].* 1/y_res*1/x_res
        integral_f_phi_1[:,(chunk*total_shape_info[2]+1):(chunk*total_shape_info[2]+data_shape[5])] = integral_f_phi_1_chunk
    end
end

if record_moments
    integral_2ndmom_1 *= length_domain
    integral_2ndmom_2 *= length_domain
end


if record_moments
    t_1_integral = mean_1[1]*height_domain*length_domain
else
    mean_1 = [sum(integral_f_phi_1[:,i].*psi_spacing_1) for i=1:nt]
    integral_2ndmom_1 = [sum(integral_f_phi_1[:,t].*psi_spacing_1.^2) for t in 1:nt]*length_domain*height_domain
    t_1_integral = sum(integral_f_phi_1[:,1].*psi_spacing_1)*length_domain*height_domain
end

plt1=plot()
plt1_1=plot()
plt1_12=plot()
plt1_2=plot()
plt2=plot()
plt3=plot()
plt4=plot()
plt5=plot()

# omega_bar = (20.0) #mean turb freq

# # u_mean = (0) #velocity mean
# turb_k_e= (5) #turbulent kinetic energy = 3/2 Var(u)
# C_0 = (5) # ~Kolmogorov constant - a rate for velocity change

# β=0.75*C_0*omega_bar
# D  = (8/9)*turb_k_e/(C_0*omega_bar)

# chi_exact_t = (mean_1[2]*length_domain).^2/(8*sqrt(2*pi*D))*t_space.^(-1.5)
# constant_wrtD(D) = integral_2ndmom_1[1]/(8*sqrt(2*pi*D))

# chi_t=-0.5*(integral_2ndmom_1[3:end]-integral_2ndmom_1[1:end-2])/(2*dt)

# max_m=argmax(chi_t)+10

# plot!(plt1_12,t_space[2:end-2],chi_t,label="χ")
# # plot!(plt1_12,t_space[max_m:end],chi_exact_t[max_m:end],label="χ theoretical")

# end_time=nt-2
# # try
# #     global end_time=minimum(findall(chi_t[max_m:end].<0))-1
# # catch e
# #     isa(e,ArgumentError) || rethrow()
# # end

# # conv_to_scientific(x) = @sprintf("%.0e", x)

# # lin_t = t_space[max_m+1:end_time+1]
# # log_t = log.(lin_t)
# # t_min =  t_space[max_m+1]
# # t_max = t_space[end_time+1]
# # t_ticks = [j*exp10(i) for j=1:9,i=floor(Int,log10(t_min)):ceil(Int,log10(t_max)) ]
# # t_ticks=round.(t_ticks,sigdigits=1)
# # in_range_t = (t_ticks.>= t_min )
# # t_ticks_locs=log.(t_ticks[in_range_t])
# # t_ticks = string.(t_ticks)
# # t_ticks[[3,5,7,9],:,:] .= ""
# # t_ticks=t_ticks[in_range_t]

# # y_max = max(maximum(chi_exact_t[max_m:end_time]),maximum(chi_t[max_m:end_time]))
# # y_min = min(minimum(chi_exact_t[max_m:end_time]),minimum(chi_t[max_m:end_time]))
# # y_ticks = [j*exp10(i) for j=1:9,i=floor(Int,log10(y_min)):ceil(Int,log10(y_max)) ]
# # y_ticks=round.(y_ticks,sigdigits=1)
# # in_range_y = (y_ticks.>= y_min ) .& (y_ticks.<=y_max )
# # y_ticks_locs=log.(y_ticks[in_range_y])
# # y_ticks = conv_to_scientific.(y_ticks)
# # y_ticks[[3,5,7,9],:,:] .= ""
# # y_ticks=y_ticks[in_range_y]

# # power, constant = linregress(log_t,log.(chi_t[max_m:end_time])).coeffs

# # println("power=",power)

# # scatter!(plt1_2,log_t,log.(chi_t[max_m:end_time]),
# # title="Evolution of scalar disipation rate over time",
# # ylabel="χ",
# # xlabel="time",
# # label="χ evaluated from simulation",
# # xticks=(t_ticks_locs,t_ticks),
# # yticks =(y_ticks_locs,y_ticks),
# # formatter = :scientific where N
# # )

# # plot!(plt1_2,log_t,log.(t_space[max_m:end_time]).*power.+constant,
# # label="Linear fit to data in log-log space")

# # plot!(plt1_2,log_t,log.(chi_exact_t[max_m:end_time]),
# # label="χ theoretical for 1D diffusion")

# # record_moments|| (mean_1=sum(mean_field_1,dims=[1,2])*1/(y_res*x_res))
# # record_moments|| (mean_2=sum(mean_field_2,dims=[1,2])*1/(y_res*x_res))
# # mass_loss_rate_1 = (mean_1[1:end-2]-mean_1[3:end])/(2*dt)
# # mass_loss_rate_2 = (mean_2[1:end-2]-mean_2[3:end])/(2*dt)
# # scatter!(plt2,t_space[2:end-2],mass_loss_rate_1, label="R1")
# # scatter!(plt2,t_space[2:end-2],mass_loss_rate_2, label="R2")
# # max_r1=argmax(mass_loss_rate_1[10:end])+9
# # plot!(plt2,t_space[max_r1+1:end],t_space[max_r1+1:end].^(-1.5).*(mass_loss_rsate_1[max_r1]./ t_space[max_r1+1].^(-1.5)))

heatmap!(plt2,t_space[1:end-1],psi_spacing_1,integral_dist_1,
levels=200, 
color= colormap("RdBu", logscale=true),
title="ϕ1 dist over time",
ylabel="ϕ",
xlabel="time",
colorbartitle="Relative frequency",
fill=true,
line=false)

# y_integral_means_1 = sum(mean_field_1, dims=1)[1,:,:] .* 1/y_res
# heatmap!(plt3,t_space[1:end-1],x_space,y_integral_means_1,
# levels=200, 
# color= colormap("RdBu", logscale=true),
# title="Evolution of mean concentration ϕ1 vs y over time",
# ylabel="x",
# xlabel="time",
# colorbartitle="concentration",
# fill=true,
# line=false)
x_integral_means_1 = sum(mean_field_1, dims=2)[:,1,:] .* 1/x_res
heatmap!(plt3,t_space[1:end-1],y_space,x_integral_means_1[:,:],
levels=200, 
color= colormap("RdBu", logscale=true),
title="Evolution of mean concentration ϕ1 vs y over time",
ylabel="y",
xlabel="time",
colorbartitle="concentration",
fill=true,
line=false)

# # y_integral_means_2 = sum(mean_field_2, dims=1)[1,:,:] .* 1/y_res
# # heatmap!(plt4,t_space[1:end-1],x_space,y_integral_means_2[:,:],
# # levels=200, 
# # # color= colormap("RdBu", logscale=true),
# # title="Evolution of mean concentration ϕ2 vs y over time",
# # ylabel="x",
# # xlabel="time",
# # colorbartitle="concentration",
# # fill=true,
# # line=false)
# heatmap!(plt4,x_space[:],y_space,mean_field_1[:,:,end],#change the final index in phi_1_means to choose time plotted
# levels=200, 
# color= colormap("RdBu", logscale=true),
# title="Mean ϕ1 vs space at the final time",
# ylabel="y",
# xlabel="x",
# colorbartitle="concentration",
# fill=true,
# line=false)

# # var_from_k=sum(prod_field, dims=[1])[1,:,:]./y_res.-0.125
# x_space_pdf=y_integral_means_1./(sum(y_integral_means_1,dims=[1]).*length_domain/x_res)
# mode_index_x = getindex.(argmax(x_space_pdf, dims=[1]),1)[:]
# #centering by mode to account for periodicity, be careful with, only works if modes are transported
# x_space_pdf_centred = x_space_pdf#zeros(size(x_space_pdf))
# # for i in 1:size(x_space_pdf)[2]
# #     x_space_pdf_centred[:,i]=circshift(x_space_pdf[:,i],(mode_index_x[1]-mode_index_x[i]))
# # end
# mean_x = sum(x_space_pdf_centred.*(x_space), dims=[1])[1,:].*length_domain/x_res
# # mean_x=length_domain/2 .+2*t_space[1:end-1]
# var_x= sum(x_space_pdf_centred.*(x_space.-transpose(mean_x)).^2, dims=[1])[1,:].*length_domain/x_res
# var_x_exact(t) = 2*D*(t-1/β*(1-exp(-β*t)))+(length_domain/50)^2
# var_x_bias=var_x[1].-(length_domain/50)^2
# var_x .-= var_x_bias
var_1 = integral_2ndmom_1/(length_domain) .- mean_1.^2
var_1_2 = [sum((psi_spacing_1.-mean_1[i]).^2 .*integral_dist_1[:,i]) for i=axes(integral_dist_1,2)]
plot!(plt4,t_space[1:end-1],var_1
,ylabel="Var(X)" 
,xlabel="t"
# ,label="Particle Solution"
#,title="Variance comparison to exact solution"
)
# plot!(plt5,t_space[1:end-1],var_x .-var_x_bias,
#     label="Particle Solution with Shifted var(X0)")
plot!(plt5,t_space[1:end-1],mean_1)#.*10/mean_1[1])
# plot!(plt5,t_space[2:end-1],-(mean_1[2:end]-mean_1[1:end-1])/(dt))
# plot!(plt5,t_space[1:end-1],var_x_exact,
#     label="Exact Solution")
# heatmap!(plt5,t_space[1:end-1],x_space,var_from_k,
#     levels=200, 
#     color= colormap("RdBu", logscale=true),
#     title="variation from constant k vs y over time",
#     ylabel="x",
#     xlabel="time",
#     colorbartitle="concentration",
#     clims=(min(-maximum(var_from_k),minimum(var_from_k)),max(maximum(var_from_k),-minimum(var_from_k))),
#     fill=true,
#     line=false)
# for y in y_space_edges
#     plot!(plt4, [0,length_domain],[y,y],label="",color="cyan")
# end
# for x in x_space_edges
#     plot!(plt4, [x,x], [0,height_domain],label="",color="cyan")
# end

# savefig(plt1,image_folder*"mean_over_time")
# savefig(plt1_1,image_folder*"rate_over_time")
# savefig(plt1_12,image_folder*"rvchi_over_time")
# savefig(plt1_2,image_folder*"chi_over_time")
# savefig(plt2,image_folder*"pdf_over_time")
# savefig(plt3,image_folder*"mean_y_over_time")
# savefig(plt4,image_folder*"initial_time_concentarions")
# savefig(plt5,image_folder*"final_time_concentarions")

# uncommentfor display. If Plots backend is set to plotlyjs these are interactive
# display(plt1) 
# display(plt1_2)
# display(plt1_12)
display(plt2)
display(plt3)
display(plt4)
display(plt5)
# #used for holding host process for plots open - not needed for vscode operation
# print(stdout, "press enter to close plots")
# read(stdin, 1)