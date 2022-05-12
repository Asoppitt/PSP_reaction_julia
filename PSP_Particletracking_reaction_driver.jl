import Random
include("PSP_Particletracking_module.jl")
pptr = PSPParticleTrackingReactions

Random.seed!(12345) #setting a seed

folders = ["Data","default","F32"]#note, would advise using file names or folders to store parameter sets, data file does not

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

T=0.10 #total time

#geometry controls
#spacial dimentions
length_domain = float_type(0.2) #length of periodic element
height_domain = float_type(0.01)
#dimesions in number of cells
n=5#a scale for resolution:use is optional, can set x_res,y_res directly
x_res=4*n #number of cells in x dim
y_res=n #number of cells in x dim

np = x_res*y_res*1000 # number of particles
#concentration space resolution for stored pdfs
psi_partions_num = 50
psi_domain = [float_type(0.0),float_type(1.01)]
#controls for omega distribution
#omega follows:
#dω=-ω(ω-̄ω)dt+√(2σ^2(ω-ω_min)ω)dW
omega_bar = float_type(10.0) #mean turb freq
omega_sigma_2 = float_type(1) #variance of turb freq - σ^2
omega_min = float_type(0.00) # minimum value for omega - acts as a shift on gamma distribution. scales as O(1/λ_max^2) where λ_max is the maximum lengthscale 
#mixing paramters
B_format = "Constant"
c_phi = float_type(25.0)#streangth of mixing
c_t = float_type(2.0)#re-pairing timescale parameter
#velocity mean and variance
u_mean = float_type(3) #velocity mean
turb_k_e= float_type(0.938) #turbulent kinetic energy = 3/2 Var(u)
C_0 = float_type(2.1) # ~Kolmogorov constant - a rate for velocity change
#strength of boundary reaction: 
#Ddϕ/dx=bc_k*\phi
bc_k = float_type(0.0)
reacting_boundaries=["lower",]
#parmeters for random bc
num_vp = float_type(Inf)
bc_CLT = true

bulk_k=float_type(0.0)
bulk_reaction=((a,b)->float_type(0),(a,b)->float_type(0))

PSP_on = true
#intial condtion one of:
#"Uniform phi_1","triple delta","2 layers","double delta","centred normal","centred 2 normal",
#("double delta difference",equlibrium_product),("2 layers difference",equlibrium_product),     these are for the bulk reaction
#("1 layer transport, 1 layer empty",empty_layer)  empty_layer=0 gives no mass at bottom, empty_layer=1 gives no mass at top
ic = ("1 layer transport, 1 layer empty",1)

nt = max(50,ceil(Int, T*3*sqrt((2/3)*turb_k_e)*max(x_res/length_domain,y_res/height_domain)))   # number of time steps
dt = float_type(T/nt)#float_type(0.0002)  # time step
println(nt)
step_95=2*sqrt((2/3)*turb_k_e)*dt#95% of steps will be less than this
(step_95>length_domain/x_res || step_95>height_domain/y_res) && @warn "cell size smaller than >5% steps"

space_cells = pptr.cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = pptr.psi_grid(psi_partions_num, psi_domain)
mix_params, move_params, bc_params = pptr.PSP_motion_bc_params(omega_bar, omega_sigma_2,omega_min, C_0, B_format, c_phi, c_t, u_mean, bc_k, num_vp, bc_CLT, bulk_reaction=bulk_reaction, reacting_boundaries=reacting_boundaries)

(bc_params.bc_k*sqrt(pi*bc_params.B/(bc_params.C_0*turb_k_e))>1) && throw(ErrorException("reaction prob >1"))

xp = zeros(float_type, np,nt+1)
yp = zeros(float_type, np,nt+1)
xp[:,1] = length_domain.*rand(float_type, np)
yp[:,1] = height_domain.*rand(float_type, np)
bc_interact = pptr.particle_motion_model(xp,yp,turb_k_e,move_params,dt,space_cells)

f_phi = zeros(float_type, psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
write(base_filename*"array_shape", [psi_partions_num, psi_partions_num, y_res, x_res, nt])
if PSP_on 
    pptr.PSP_model!(f_phi ,xp, yp, turb_k_e, bc_interact, dt, ic, mix_params, psi_mesh, space_cells, bc_params, verbose)
else
    pptr.make_f_phi_no_PSP!(f_phi,xp,yp, turb_k_e, bc_interact, ic, psi_mesh, space_cells, bc_params, verbose)
end
write(base_filename*"data", f_phi)
println("data saved at:"*base_filename*"data")