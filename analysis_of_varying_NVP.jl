using Plots; plotlyjs()
base_filename = "Data/PSP_on_uniform_1_square_c_0_21_k_09_w_1_new_abs"

n = 10
nt = 60
psi_partions_num = 20
phi_domain = [0,1.2]
Del_phi = (phi_domain[2]-phi_domain[1])/psi_partions_num
Del_phi_x_phi = (Del_phi)^2 
x_res=n #number of cells in x dim
y_res=n
in_data=zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)

binomial_NVP = Array(1:20)
append!(binomial_NVP,Array(22:2:78))
append!(binomial_NVP, Array(2 .^(3:10).*10))

CLT_NVP = Array(1:20)
append!(CLT_NVP,Array(22:2:78))
append!(CLT_NVP, Array(2 .^(3:10).*10))
# binomial_NVP = [1:20;22:2:40]
# CLT_NVP = binomial_NVP
# binomial_NVP = 1:4
# CLT_NVP = binomial_NVP

psi_spacing = LinRange(0,1.2,psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

final_time_xedge_bino = zeros(psi_partions_num,size(binomial_NVP)[1])

for (NVP,i) in zip(binomial_NVP,1:size(binomial_NVP)[1])
    filename = base_filename*"_"*string(NVP)*"_vp"
    read!(filename,in_data)
    final_time_xedge_bino[:,i] = sum(in_data,dims=[2,4])[:,1,1,1,nt].* 1/x_res./Del_phi
end 

final_time_xedge_CLT = zeros(psi_partions_num,size(CLT_NVP)[1])

for (NVP,i) in zip(CLT_NVP,1:size(CLT_NVP)[1])
    filename = base_filename*"_"*string(NVP)*"_vp_CLT"
    read!(filename,in_data)
    final_time_xedge_CLT[:,i] = sum(in_data,dims=[2,4])[:,1,1,1,nt].* 1/x_res./Del_phi#Del_phi normalises to an actual pdf
end 
final_time_xedge_inf = zeros(psi_partions_num)
for i in [1]
    filename = base_filename*"_Inf_vp_CLT"
    read!(filename,in_data)
    global final_time_xedge_inf = sum(in_data,dims=[2,4])[:,1,1,1,nt].* 1/x_res./Del_phi
end
conv_error_clt_to_mean=zeros(size(CLT_NVP)[1])
for i in 1:size(CLT_NVP)[1]
    conv_error_clt_to_mean[i] = sum(abs.(final_time_xedge_CLT[:,i]-final_time_xedge_inf).^2)*Del_phi
end 
conv_error_bino_to_mean=zeros(size(binomial_NVP)[1])
for i in 1:size(binomial_NVP)[1]
    conv_error_bino_to_mean[i] = sum(abs.(final_time_xedge_bino[:,i]-final_time_xedge_inf).^2)*Del_phi
end 
n_clt_and_bino = argmax(CLT_NVP[in.(CLT_NVP,[binomial_NVP])])
conv_error_bino_to_clt=zeros(n_clt_and_bino)
for i in 1:n_clt_and_bino
    conv_error_bino_to_clt[i] = sum(abs.(final_time_xedge_CLT[:,i]-final_time_xedge_bino[:,binomial_NVP.==CLT_NVP[i]]).^2)*Del_phi
end 

plt1 = scatter(log.(CLT_NVP),conv_error_clt_to_mean)
scatter!(log.(binomial_NVP[1:end]),conv_error_bino_to_mean[1:end])
ticks_locs = 2 .^(floor(log2(minimum(CLT_NVP))):log2(maximum(CLT_NVP)))
xticks!(log.(ticks_locs),string.(ticks_locs))
xlabel!(plt1,"Number of Virtual Particles")
ylabel!(plt1,"Error in L2 norm")
title!(plt1,"Error of Binomial and Normal compared to Mean model")
display(plt1)
plt2 = plot(scatter(log.(CLT_NVP[1:n_clt_and_bino]),conv_error_bino_to_clt))
N_range = Array(CLT_NVP[1]:CLT_NVP[n_clt_and_bino])
plot!(log.(N_range),1 ./sqrt.(2 .*pi.*log.(N_range)))
xlabel!(plt2,"Number of Virtual Particles")
ylabel!(plt2,"Error in L2 norm")
title!(plt2,"Error of Binomial compared to Normal")
display(plt2)
plt3=surface(log.(binomial_NVP),psi_spacing,final_time_xedge_bino, color= colormap("RdBu", logscale=true))
display(plt3)
surface(log.(CLT_NVP),psi_spacing,final_time_xedge_CLT, color= colormap("RdBu", logscale=true))# xaxis=log)