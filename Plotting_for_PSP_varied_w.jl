using Plots; GR.init()#plotlyjs()
import Statistics as st

plt=plot()
base_filename = "Data/test_lang_var_w_edge_data/F32/"
data_shape = zeros(Int,5)
read!(base_filename*"array_shape",data_shape)
#  [psi_partions_num, psi_partions_num, y_res, x_res, nt+1]
# data_shape = [50,50,5,15,114]
nt = data_shape[5]
psi_partions_num = data_shape[1]
phi_domain = [0,1.01]
x_res=data_shape[4] #number of cells in x dim
y_res=data_shape[3]
plt=plot()
plt_1=plot()
plt_2=plot()
for omega_bar =  Float32.([4,12,16,22,24,28,32,40,50])
    for NVP= [Inf]
        (NVP<Inf) && (NVP=Int(NVP))
        filename=base_filename*string(NVP)*"_vp"
        (NVP==1)||(filename*="_CLT")
        filename*="_w"*string(omega_bar)
        in_data=zeros(Float32,psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
        # in_data_flux=zeros(2,x_res, nt)
        # in_data_flux_vel=zeros(x_res, nt)
        length_domain = 0.03 #length of periodic element
        height_domain = 0.01
        T=0.6
        dt=0.2/(nt)
        k = 0.938
        C_0 = 2.1
        omega_mean=10
        B=(1.5*C_0)
        bc_k=0.25
        D=2*C_0*k/((B^2)*omega_mean)
        read!(filename,in_data)
        edge_mean=zeros(Float32,nt+1)
        edge_2=zeros(Float32,nt+1)
        edge_2_v=zeros(Float32,nt)
        read!(filename*"edge_mean", edge_mean)
        read!(filename*"edge_squared", edge_2)
        read!(filename*"edge_squared_velocity", edge_2_v)
        t_space = LinRange(0,nt*dt,nt+1)
        x_space = LinRange(0,length_domain,x_res)
        y_space = LinRange(0,height_domain,y_res)

        psi_spacing = LinRange(phi_domain[1],phi_domain[2],psi_partions_num+1) #defines boundaries for cells in psi space
        psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

        phi_1_means = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]
        phi_2_means = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing ) for psi_1_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]

        x_integral_means = sum(phi_1_means, dims=2)[:,1,:] .* 1/x_res*length_domain
        integral_mean_1 = sum(x_integral_means, dims=1)[1,:] .* 1/y_res*height_domain
        y_integral_means = sum(phi_1_means, dims=1)[1,:,:] .* 1/y_res*height_domain

        x_integral_f_phi_1 = sum(in_data,dims=[2,4])[:,1,:,1,:].* 1/x_res#*length_domain
        integral_f_phi_1 = sum(x_integral_f_phi_1,dims=2)[:,1,:].* 1/y_res#*height_domain
        integral_phi_1_2nd_mom = [sum(integral_f_phi_1[:,t].*psi_spacing.^2 ) for t in 1:nt+1]*length_domain*height_domain
        edge_2*=length_domain
        edge_2_v*=length_domain

        plt1=plot()
        plt1_1=plot()
        plt1_12=plot()
        plt1_2=plot()
        plt1_3=plot()
        plt2=plot()
        
        plot!(plt1,t_space,integral_mean_1,
        title="Evolution of mean concentration over time",
        ylabel="concentration",
        xlabel="time",
        label=""
        )

        dmdt=-0.5*(integral_mean_1[3:end]-integral_mean_1[1:end-2])/dt
        dc_2dt = -0.5*(integral_phi_1_2nd_mom[3:end]-integral_phi_1_2nd_mom[1:end-2])/dt

        plot!(plt1_1,t_space[2:end-1],dc_2dt,
        title="dc^2/dt reaction rate over time",
        ylabel="reaction rate",
        xlabel="time",
        label=""
        )

        # plot!(plt1_12,t_space[2:end-1],dmdt)
        plot!(plt1_12,t_space[1:end],edge_2)
        plot!(plt1_12,t_space[2:end-1],dc_2dt)
        plot!(plt1_2,t_space[1:end],edge_2,label="e^2")

        G=2*omega_bar/2

        # plot!(plt1_3,t_space[2:end-1], -0.5*(dc_2dt-edge_2[2:end-1]))
        plot!(plt,(t_space[2:end-1]), (-0.5*(dc_2dt-G*edge_2[2:end-1])), label=omega_bar)
        plot!(plt_1,(t_space[2:end-1]), ((dc_2dt)), label=omega_bar)
        plot!(plt_2,(t_space[2:end-1]), (G*edge_2[2:end-1]), label=omega_bar)
        # plot!(plt1_3,t_space[2:end-1], -0.5*(dc_2dt-edge_2[2:end-1]+edge_2_v[2:end]))

        # plot!(plt1_12,t_space[2:end-1],dmdt,
        # title="Comparison of scalar disipation rate and reaction rate",
        # ylabel="reaction rate",
        # xlabel="time",
        # label=""
        # )
        # plot!(plt1_12,t_space[2:end-1],-0.5*(integral_2ndmom_1[3:end]-integral_2ndmom_1[1:end-2])/dt,label="")

        # plot!(plt1_2,t_space[2:end-1],-0.5*(integral_2ndmom_1[3:end]-integral_2ndmom_1[1:end-2])/dt,
        # title="Evolution of scalar disipation rate over time",
        # ylabel="χ",
        # xlabel="time",
        # label=""
        # )

        heatmap!(plt2,t_space[:],psi_spacing[1:ceil(Int,end)],integral_f_phi_1[1:ceil(Int,end),:], levels=100,
        color= colormap("RdBu", logscale=true),
        # title="Contour Plot of ϕ denisty over time for N → ∞",
        ylabel="Normalised Concentration",
        xlabel="Time (dimensionless)",
        # guidefontsizes=18,
        # colorbar_titlefontsize=18,
        # tickfontsize=14,
        # clims=[0,12],
        # colorbar=:none,
        colorbartitle="Relative frequency" )
        # contour!(plt2,t_space[:],psi_spacing,integral_f_phi_1[:,:],
        # levels=200, 
        # color= colormap("RdBu", logscale=true),
        # # title=string(NVP),
        # ylabel="ϕ",
        # xlabel="time",
        # colorbartitle="Relative frequency",
        # fill=true,
        # line=false)
        # plt=plot(psi_spacing,integral_f_phi_1[:,Int(floor(end/4))], ylabel="Relative frequency", xlabel="ϕ",
        # # title="pdf of ϕ at time 0.3 for N = 1 with PSP mixing",
        # label="")
        # savefig("dist_plots/thin"*string(NVP)*"quartertime")
        # plot(t_space,x_integral_means[1,:])
        # display(plt1)
        # display(plt1_1)
        # display(plt1_12)
        # display(plt1_2)
        # display(plt1_3)
        # display(plt2)
        try
            folder = "poster_plots/var_w/"
            run(`mkdir $folder`)
        catch e 
            isa(e,ProcessFailedException) || rethrow()
        end
        # savefig(plt1,"poster_plots/var_w/"*string(Int(omega_bar)))#*string(NVP)*"NVP_")
        # savefig(plt2,"dist_plots/thin"*string(NVP)*"contour_betterdecorr")
    end
end
display(plt)
display(plt_1)
display(plt_2)
# savefig(plt, "dist_plots/quartertime_stack")