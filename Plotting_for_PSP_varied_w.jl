using Plots; GR.init()#plotlyjs()
import Statistics as st

plt=plot()
base_filename = "Data/poster/var_NVP/F32/"
data_shape = zeros(Int,5)
read!(base_filename*"array_shape",data_shape)
#  [psi_partions_num, psi_partions_num, y_res, x_res, nt+1]
# data_shape = [50,50,5,15,114]
nt = data_shape[5]
psi_partions_num = data_shape[1]
phi_domain = [0,1.01]
x_res=data_shape[4] #number of cells in x dim
y_res=data_shape[3]
for omega_bar =  Float32.([8])#[4,8,12,16,18,20])#[1,4,8,12,16,18,20,22,24,28,32,40,50])
    for NVP= [50,200,Float32(Inf)]
        (NVP<Inf) && (NVP=Int(NVP))
        filename=base_filename*string(NVP)*"_vp"
        (NVP==1)||(filename*="_CLT")
        filename*="_w"*string(omega_bar)
        in_data=zeros(Float32,psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
        # in_data_flux=zeros(2,x_res, nt)
        # in_data_flux_vel=zeros(x_res, nt)
        length_domain = 0.20 #length of periodic element
        height_domain = 0.05
        T=0.6
        dt=0.2/(nt)
        k = 0.938
        C_0 = 2.1
        omega_mean=5
        B=(1.5*C_0)
        bc_k=0.25
        D=2*C_0*k/((B^2)*omega_mean)
        do_flux = false
        read!(filename,in_data)
        do_flux && read!(filename*"flux",in_data_flux)
        do_flux && read!(filename*"flux_vel",in_data_flux_vel)
        t_space = LinRange(0,nt*dt,nt+1)
        x_space = LinRange(0,length_domain,x_res)
        y_space = LinRange(0,height_domain,y_res)

        psi_spacing = LinRange(phi_domain[1],phi_domain[2],psi_partions_num+1) #defines boundaries for cells in psi space
        psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])#defining centres of bins

        phi_1_means = [sum([sum(in_data[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]
        phi_2_means = [sum([sum(in_data[psi_1_index,:,i,j,t].*psi_spacing ) for psi_1_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt+1]

        x_integral_means = sum(phi_1_means, dims=2)[:,1,:] .* 1/x_res
        y_integral_means = sum(phi_1_means, dims=1)[1,:,:] .* 1/y_res

        x_integral_f_phi_1 = sum(in_data,dims=[2,4])[:,1,:,1,:].* 1/x_res
        integral_f_phi_1 = sum(x_integral_f_phi_1,dims=2)[:,1,:].* 1/y_res

        plt1=plot()
        plt2=plot()
        # do_flux && plot!(t_space,flux_y0_1_int)
        # do_flux && plot(t_space,flux_y0_vel)
        # surface(t_space,y_space,y_integral_means,ylabel="y",
        # xlabel="t",title="Particle method")
        # plot!(plt,t_space,[ x_integral_means[1,i] for i=1:nt+1])
        # plot!(t_space,[ x_integral_means[end,i] for i=1:nt+1])
        # plot!(t_space,[ st.mean(x_integral_means[:,i]) for i=1:nt+1])
        # surface(x_space,y_space,phi_1_means[:,:,nt], ylabel="y")
        contourf!(plt1,t_space[:],psi_spacing[1:ceil(Int,end)],integral_f_phi_1[1:ceil(Int,end),:], levels=100,
        color= colormap("RdBu", logscale=true),
        # title="Contour Plot of ϕ denisty over time for N → ∞",
        ylabel="ϕ",
        xlabel="time",
        colorbartitle="Relative frequency")
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
        display(plt1)
        # display(plt2)
        # try
        #     folder = "poster_plots/var_NVP/"
        #     run(`mkdir $folder`)
        # catch e 
        #     isa(e,ProcessFailedException) || rethrow()
        # end
        savefig(plt1,"poster_plots/var_NVP/"*string(NVP)*"contour")
        # savefig(plt2,"dist_plots/thin"*string(NVP)*"contour_betterdecorr")
    end
end
# savefig(plt, "dist_plots/quartertime_stack")