import Random
using TurbulentMixingParticleTrackingReactions
using ProfileView
using Distributions
# include("/home/u2093579/Documents/Project_code/TurbulentMixingParticleTrackingReactions/src/PSP_Particletracking_module.jl")
# include("/home/u2093579/Documents/Project_code/TurbulentMixingParticleTrackingReactions/src/save_as_you_go.jl")

Random.seed!(12345) #setting a seed

folders = ["Data","test_0","1000_part","F32"]#note, would advise using file names or folders to store parameter sets, data file does not

base_filename=""
for folder in folders
    global base_filename*=folder
    try
        run(`mkdir $base_filename`)
    catch e 
        isa(e,ProcessFailedException) || rethrow()
    end
    global base_filename*="/"
end

float_type = Float32

verbose=true #printing of step numbers

T=250000*7.63364e-09 #total time - is not stored as part of data
nt=2
#geometry controls
#spacial dimentions - only rectangular geometry is simulated
length_domain = float_type(5e-3) #length of periodic element
@show height_domain = float_type(length_domain.*1/8)#float_type(0.05)
#dimesions in number of cells
n=20#a scale for resolution:use is optional, can set x_res,y_res directly
x_res=8*n #number of cells in x dim
y_res=n #number of cells in  dim

dim_phi=float_type(10/10)
dim_t=float_type(T /T)
T=T.*dim_t
length_domain_nd=length_domain#float_type(2)
dim_x=float_type(length_domain_nd/length_domain)#going to try non-dimensionalising space
length_domain*=dim_x
# # n=2
height_domain = float_type(length_domain.*0.5)

np = x_res*y_res*1000 # number of particles
#concentration space resolution for stored pdfs
psi_partions_num_1 = min(100,floor(Int,sqrt(np//x_res*y_res)))
psi_partions_num_2 = 2# min(100,floor(Int,sqrt(np//x_res*y_res)))   
psi_domain = [float_type(0.0),float_type(1.01)]*dim_phi
#controls for omega distribution
#omega follows:
#dω=-ω(ω-̄ω)dt+√(2σ^2(ω-ω_min)ω)dW
#giving statistics E[ω]=omega_bar, Var[ω]=omega_sigma_2
a=LogNormal(2.75384,2.0475)#(3.8975985750068634,2.8665512122893175)#(3.8731577179478354,2.585404654956942)
omega_bar = float_type(mean(a))./dim_t #mean turb frequency
omega_sigma_2 = float_type(var(a))./dim_t^2#float_type(1.60090716183498e11)# #variance of turb freq - σ^2
omega_min = float_type(0.00)./dim_t # minimum value for omega - acts as a shift on gamma distribution. scales as O(1/λ_max^2) where λ_max is the maximum lengthscale 
omega_dist=:LogNormal#:Gamma#:LogNormal
a=nothing
#mixing paramters
c_phi = float_type(25.0)#streangth of mixing
c_t = float_type(2.0).*dim_t#re-pairing timescale parameter
#velocity mean and variance
#velocity follows:
#u=u_mean+u^*
#du^*=(-0.5*B(C_0)*omega_bar*u^*)*dt+sqrt.(C_0.*turb_k_e[:,t].*omega_bar)dW
#giving statistics E[u]=u_mean, Var[u]=turb_k_e
u_mean = float_type(0.36)*dim_x/dim_t #velocity mean

if false # reading mean velocity from file
    solid_thickness=max(ceil(Int,3000/2),4)
    xdata_in=read("Data/DNS/velocity_x5000.dat", String)
    # xdata_in=xdata_in[1:end-2]#strip trailing space and newline
    xdata = split(xdata_in, ' ')
    ux_data = parse.(float_type,reshape(xdata[xdata .!= "\n"],(3000,3000)))
    ux_data = ux_data[solid_thickness:end,:]
    ydata_in=read("Data/DNS/velocity_y5000.dat", String)
    # ydata_in=ydata_in[1:end-2]#strip trailing space and newline
    ydata = split(ydata_in, ' ',)
    uy_data =  reshape(parse.(float_type,ydata[ydata .!= "\n"]),(3000,3000))
    uy_data = uy_data[solid_thickness:end,:]
    u_mean = ((x,y)-> ux_data[floor(Int,x/length_domain*(size(ux_data,1)-1))+1,floor(Int,y/height_domain*(size(ux_data,2)-1)+1)],
            (x,y)-> uy_data[floor(Int,x/length_domain*(size(uy_data,1)-1))+1,floor(Int,y/height_domain*(size(uy_data,2)-1)+1)])
end

turb_k_e= float_type(0.18785770884020575)*(dim_x/dim_t)^2 #turbulent kinetic energy = 3/2 Var(u)
C_0 = float_type(5) # ~Kolmogorov constant - a rate for velocity change
B_format = "Constant"
#strength of boundary reaction: 
#Ddϕ/dx=bc_k*ϕ at boundary
@show bc_k =  float_type(2190)#float_type(2)*dim_x/(dim_t)#float_type(79418.40738854579) #strength of boundary reaction
reacting_boundaries=["lower",]#upper,lower,left,right
#parmeters for random bc
num_vp = 20#float_type(Inf)
lang_half_sat=float_type(0.75)#langmuir saturation constant, value of phi for which removal rate is halfed
# bc_non_lin_corr(phi) = lang_half_sat ./(lang_half_sat .+phi) #langmuir type correction
bc_non_lin_corr=nothing
# bc_non_lin_corr(phi) = [phi[i]==0 ? 0 : 1 ./sqrt.(phi[i]) for i=1:2]

bulk_k=float_type(0.0)
rate=float_type(100)
# bulk_reaction=((a,b)->(0.5*(a-b+sqrt((a-b)^2+4*bulk_k))-a)/dt,(a,b)->(0.5*(b-a+sqrt((a-b)^2+4*bulk_k))-b)/dt) #"instant" reaction
# bulk_reaction=((a,b)->-rate*(bulk_k+a*b),(a,b)->-rate*(bulk_k+a*b))
bulk_reaction=((a,b)->-0,(a,b)->-0)

PSP_on = true
#intial condtion one of:
#"Uniform phi_1","triple delta","2 layers","double delta","centred normal","centred 2 normal","empty"
#("double delta difference",equlibrium_product),("2 layers difference",equlibrium_product),("vertical strip difference",left_edge,right_edge,equlibrium_product) these are for the bulk reaction
#("vertical strip",left_edge,right_edge)
#("1 layer transport, 1 layer empty",empty_layer)  empty_layer=0 gives no mass at bottom, empty_layer=1 gives no mass at top
ic = ("Uniform phi_1",float_type(1 .*dim_phi))#("x centred normal",length_domain/float_type(50))

nt = max(nt,ceil(Int, T*3*sqrt((2/3)*turb_k_e)*max(x_res/length_domain,y_res/height_domain)))   # number of time steps
dt = float_type(T/nt)#float_type(0.0002)  # time step
println(nt)
step_95=2*sqrt((2/3)*turb_k_e)*dt#95% of steps will be less than this
(step_95>length_domain/x_res || step_95>height_domain/y_res) && @warn "cell size smaller than >5% steps"

space_cells = cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = psi_grid(psi_partions_num_1,psi_partions_num_2, psi_domain)
mix_params, move_params, bc_params = PSP_motion_bc_params(omega_bar, omega_sigma_2,omega_min, C_0, B_format, c_phi, c_t, u_mean, bc_k, num_vp, bulk_reaction=bulk_reaction, reacting_boundaries=reacting_boundaries, corr_func=bc_non_lin_corr, omega_dist=omega_dist)#,omega_t=float_type(0.1))

# no_psp_motion_model!(base_filename[1:end-1],turb_k_e, nt, dt, np, ic, move_params, mix_params, psi_mesh, space_cells, bc_params, true, 5, record_moments=true)

(bc_params.bc_k*sqrt(pi*bc_params.B/(bc_params.C_0*turb_k_e))>1) && (@warn "reaction prob >1, setting maximum realised bc_k"; @show bc_k=1/sqrt(pi*bc_params.B/(bc_params.C_0*turb_k_e)))
# error("test")
# @warn "full bc_k not realised"

if PSP_on 
    ProfileView.@profview PSP_model!(base_filename[1:end-1],turb_k_e, nt, dt, np, ic, move_params, mix_params, psi_mesh, space_cells, bc_params, true, 5; record_moments=true, saveing_rate=1) 
else
    no_psp_motion_model!(base_filename[1:end-1],turb_k_e, nt, dt, np, ic, move_params, mix_params, psi_mesh, space_cells, bc_params, true, 5, record_moments=true)
end