begin
    using Plots; GR.init()#plotlyjs()
    import Statistics as st

    plt=plot()
    plt_2nd=plot()
    base_filename_react = "Data/poster/100_phi_varw/F32/"
    base_filename_control = "Data/poster/100_phi_varw_no_border/F32/"
    data_shape = zeros(Int,5)
    read!(base_filename_react*"array_shape",data_shape)
    #  [psi_partions_num, psi_partions_num, y_res, x_res, nt+1]
    # data_shape = [50,50,5,15,114]
    nt = data_shape[5]
    psi_partions_num = data_shape[1]
    phi_domain = [0,2.01]
    x_res=data_shape[4] #number of cells in x dim
    y_res=data_shape[3]
    for omega_bar =  Float32.([4,8,12,16,18,20])#[1,4,8,12,16,18,20,22,24,28,32,40,50])
        for NVP= [200]#[50,200,Float32(Inf)]
            (NVP<Inf) && (NVP=Int(NVP))
            filename_r=base_filename_react*string(NVP)*"_vp"
            filename_c=base_filename_control*string(NVP)*"_vp"
            (NVP==1)||(filename_r*="_CLT")
            (NVP==1)||(filename_c*="_CLT")
            filename_r*="_w"*string(omega_bar)
            filename_c*="_w"*string(omega_bar)
            in_data_r=zeros(Float32,psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
            in_data_c=zeros(Float32,psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
            # in_data_flux=zeros(2,x_res, nt)
            # in_data_flux_vel=zeros(x_res, nt)
            length_domain = 0.20 #length of periodic element
            height_domain = 0.05
            T=0.2
            dt=T/(nt)
            k = 0.938
            C_0 = 2.1
            omega_mean=5
            B=(1.5*C_0)
            bc_k=0.25
            D=2*C_0*k/((B^2)*omega_mean)
            do_flux = false
            read!(filename_r,in_data_r)
            read!(filename_c,in_data_c)
            t_space = LinRange(0,nt*dt,nt+1)
            x_space = LinRange(0,length_domain,x_res)
            y_space = LinRange(0,height_domain,y_res)

            psi_spacing = LinRange(phi_domain[1],phi_domain[2],psi_partions_num+1) #defines boundaries for cells in psi space
            psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins
            

            phi_1_means_r = [sum([sum(in_data_r[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]
            phi_1_2ndmom_c = [sum([sum(in_data_c[:,psi_2_index,i,j,t].*psi_spacing.^2 ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]

            react_rates = -0.5*(phi_1_means_r[:,:,3:end]-phi_1_means_r[:,:,1:end-2])/dt
            react_rate = sum(react_rates,dims=[1,2])[1,1,:].* 1/(y_res*x_res)
            react_rate_2nd_mom = sum(react_rates.^2,dims=[1,2])[1,1,:].* 1/(y_res*x_res)
            integral_2ndmom_1_c = sum(phi_1_2ndmom_c,dims=[1,2])[1,1,:].* 1/(y_res*x_res)
            chi = -0.5*(integral_2ndmom_1_c[3:end]-integral_2ndmom_1_c[1:end-2])/dt
            chi_max=argmax(chi)
            # plt1=plot()
            plot!(plt,t_space[2:end-1],react_rate_2nd_mom-react_rate.^2)
            # plot!(plt1,t_space[2:end-1],chi)
            # display(plt1)

            # scatter!(plt,chi[chi_max:end],react_rate[chi_max:end], label=omega_bar,markeralpha=0.3,markerstrokewidth=0)
            # scatter!(plt_2nd,chi[react_rate_2nd_mom.<100],react_rate_2nd_mom[react_rate_2nd_mom.<100], label=omega_bar,markeralpha=0.3,markerstrokewidth=0)
        end
    end
    display(plt)
    # display(plt_2nd)
    # savefig(plt, "dist_plots/quartertime_stack")
end