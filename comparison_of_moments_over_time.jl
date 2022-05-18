using Plots; plotlyjs()#GR.init()
import Statistics as st

freqs=Float32.([20,50,1])#1,4,8,12,16,18,20,22,24,28,32,40,50])#[1,4]
sort!(freqs)

length_domain = 0.03#0.20 #length of periodic element
height_domain = 0.01#0.01
n=5
nt = 100#171#114#341#1137#228#682#273#137#
T=0.01
dt=T/nt
psi_partions_num = 50
phi_domain = [0,2.01]
Del_phi = (phi_domain[2]-phi_domain[1])/psi_partions_num
Del_phi_x_phi = (Del_phi)^2 
x_res=n*3#n*20 #number of cells in x dim
y_res=n#n
t_space = LinRange(0,nt*dt,nt+1)
x_space = LinRange(0,length_domain,x_res)
y_space = LinRange(0,height_domain,y_res)
psi_spacing = LinRange(phi_domain[1],phi_domain[2],psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins


base_filename = "Data/paper/toy_prob_11/with_phi_gamma/F32/"
for NVP= [200]#Float32(Inf)]
    ave_mass_var_1=zeros(Float32,length(freqs),nt+1)
    ave_mass_var_2=zeros(Float32,length(freqs),nt+1)
    ave_mass_1=zeros(Float32,length(freqs),nt+1)
    ave_mass_2=zeros(Float32,length(freqs),nt+1)
    total_rate = zeros(Float32,length(freqs))
    dphi1_dt = zeros(Float32,length(freqs),nt-1)
    dphi2_dt = zeros(Float32,length(freqs),nt-1)
    # ave_adj_x2=zeros(Float32,length(freqs),nt+1)
    plt=plot()
    plt1=plot()
    plt2=plot()
    plt3=plot()
    for (i,omega_bar) = enumerate(freqs)
        (NVP<Inf) && (NVP=Int(NVP))
        filename=base_filename*string(NVP)*"_vp"
        (NVP==1)||(filename*="_CLT")
        filename*="_w"*string(omega_bar)
        in_data=zeros(Float32,psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
        read!(filename,in_data)

        
        integral_f_phi = sum(in_data,dims=[3,4])[:,:,1,1,:].* 1/x_res* 1/y_res
        c1c2 = [sum([sum(integral_f_phi[psi_i,:,t].*psi_spacing*psi_spacing[psi_i]) for psi_i in 1:psi_partions_num])
         for t=1:nt+1]
        # x_integral_f_phi_2 = sum(in_data,dims=[1,4])[1,:,:,1,:].* 1/x_res
        integral_f_phi_2 = sum(integral_f_phi,dims=1)[1,:,:]
        # x_integral_f_phi_1 = sum(in_data,dims=[2,4])[:,1,:,1,:].* 1/x_res
        integral_f_phi_1 = sum(integral_f_phi,dims=2)[:,1,:]
        ave_mass_1[i,:]=[sum(integral_f_phi_1[:,t].*psi_spacing) for t=1:nt+1]
        ave_mass_2[i,:]=[sum(integral_f_phi_2[:,t].*psi_spacing) for t=1:nt+1]
        # ave_adj_x2[i,:]=[sum(x_integral_f_phi_1[:,1,t].*(psi_spacing).^2) for t=1:nt+1]
        ave_mass_var_2[i,:]=[sum(integral_f_phi_2[:,t].*(psi_spacing.-ave_mass_1[i,t]).^2) for t=1:nt+1]
        ave_mass_var_1[i,:]=[sum(integral_f_phi_1[:,t].*(psi_spacing.-ave_mass_2[i,t]).^2) for t=1:nt+1]
        # plot!(plt,t_space,ave_mass[i,1] .-ave_mass[i,:],label=string(omega_bar))
        plot!(plt,t_space,ave_mass_1[i,:],label="C_1 "*string(omega_bar),color=i)
        plot!(plt,t_space,ave_mass_2[i,:],label="C_2 "*string(omega_bar),color=i, linestyle=:dash)
        plot!(plt1,t_space,ave_mass_var_1[i,:],label="C_1 "*string(omega_bar),color=i)
        plot!(plt1,t_space,ave_mass_var_2[i,:],label="C_2 "*string(omega_bar),color=i, linestyle=:dash)
        dphi1_dt[i,:] = ave_mass_1[i,3:nt+1]-ave_mass_1[i,1:nt-1]
        dphi2_dt[i,:] = ave_mass_2[i,3:nt+1]-ave_mass_2[i,1:nt-1]
        plot!(plt2,t_space[2:end-1],dphi1_dt[i,:],label="C_1 "*string(omega_bar),color=i)
        plot!(plt2,t_space[2:end-1],dphi1_dt[i,:],label="C_1 "*string(omega_bar),color=i, linestyle=:dash)

        plot!(plt3,t_space,c1c2,label="C_1C_2 "*string(omega_bar),color=i)
        total_rate[i] = (ave_mass_1[i,1]-ave_mass_1[i,end])*dt*nt
        println(omega_bar,' ',ave_mass_var_1[i,1],' ', ave_mass_var_1[i,end])
    end
    # for (i,omega_bar) = enumerate(freqs)
    #     i+=length(freqs)
    #     (NVP<Inf) && (NVP=Int(NVP))
    #     filename=base_filename*string(NVP)*"_vp"
    #     (NVP==1)||(filename*="_CLT")
    #     filename*="_w"*string(omega_bar)*"_no_PSP"
    #     in_data=zeros(Float32,psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
    #     read!(filename,in_data)

    #     psi_spacing = LinRange(phi_domain[1],phi_domain[2],psi_partions_num+1) #defines boundaries for cells in psi space
    #     psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

    #     phi_1_means = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]
    #     phi_2_means = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing ) for psi_1_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]


    #     x_integral_means = sum(phi_1_means, dims=2)[:,1,:] .* 1/x_res
    #     y_integral_means = sum(phi_1_means, dims=1)[1,:,:] .* 1/y_res

    #     x_integral_f_phi_1 = sum(in_data,dims=[2,4])[:,1,:,1,:].* 1/x_res
    #     integral_f_phi_1 = sum(x_integral_f_phi_1,dims=2)[:,1,:].* 1/y_res
    #     ave_mass[i,:]=[sum(integral_f_phi_1[:,t].*psi_spacing) for t=1:nt+1]
    #     ave_mass_var[i,:]=[sum(integral_f_phi_1[:,t].*(psi_spacing.-ave_mass[i,t]).^2) for t=1:nt+1]
    # end

    # # d_ln_mean_dt = 0.5*(log.(ave_mass[:,3:end]) - log.(ave_mass[:,1:end-2]))
    # M=(ave_mass_var+ave_mass.^2)
    # # d_v_dt = 0.5*(ave_mass_var[:,3:end] - ave_mass_var[:,1:end-2])
    # d_M_dt = 0.5*(M[:,3:end] - M[:,1:end-2])
    # plt=plot(t_space[2:end-1],permutedims(0.3*ave_adj_x2[:,2:end-1]-0.5 .*d_M_dt),label="")
    # display(plt)
    # plt=plot(t_space,permutedims(M),label="")
    # display(plt)
    # plt=plot(t_space[:],permutedims(ave_adj_x2),label="")
    # display(plt)
    # plt=plot(t_space,permutedims(1 .-ave_mass),label="")
    display(plt)
    display(plt1)
    display(plt2)
    display(plt3)
    pltend = plot(freqs,total_rate)
    display(pltend)
    # ave_mass_loss=[(ave_mass[i,1] -ave_mass[i,j]) for i=1:length(freqs), j=1:nt+1]
    # ave_mass_loss_diff = [(ave_mass_loss[end,j] -ave_mass_loss[i,j]) for j=1:nt+1, i=1:length(freqs)]
    # plt_2=plot(t_space,ave_mass_loss_diff)
    # display(plt_2)
    # plt=plot(freqs,maximum(ave_mass_var, dims=2)[:,1],label="")
end