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

T=0.10
length_domain = float_type(0.2) #length of periodic element
height_domain = float_type(0.05)
n=5
x_res=4*n #number of cells in x dim
y_res=n
np = x_res*y_res*1000 # number of particles
psi_partions_num = 50
psi_domain = [float_type(0.0),float_type(1.01)]
omega_bar = float_type(10.0)
omega_sigma_2 = float_type(0.25)
C_0 = float_type(2.1)
B_format = "Constant"
c_phi = float_type(25.0)
c_t = float_type(2.0)
u_mean = float_type(0.1)
turb_k_e= float_type(1.938)
bc_k = float_type(0.3)
num_vp = float_type(Inf)
bc_CLT = true
PSP_on = true

nt = max(50,ceil(Int, T*2*sqrt((2/3)*turb_k_e)*min(x_res/length_domain,y_res/height_domain)))   # number of time steps
dt = float_type(T/nt)#float_type(0.0002)  # time step
println(nt)
step_95=2*sqrt((2/3)*turb_k_e)*dt#95% of steps will be less than this
(step_95>length_domain/x_res || step_95>height_domain/y_res) && @warn "cell size smaller than >5% steps"

space_cells = pptr.cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = pptr.psi_grid(psi_partions_num, psi_domain)
mix_params, move_params, bc_params = pptr.PSP_motion_bc_params(omega_bar, omega_sigma_2, C_0, B_format, c_phi, c_t, u_mean, bc_k, num_vp, bc_CLT)

(bc_params.bc_k*sqrt(pi*bc_params.B/(bc_params.C_0*turb_k_e))>1) && throw(ErrorException("reaction prob >1"))

xp = zeros(float_type, np,nt+1)
yp = zeros(float_type, np,nt+1)
xp[:,1] = length_domain.*rand(float_type, np)
yp[:,1] = height_domain.*rand(float_type, np)
bc_interact = pptr.particle_motion_model(xp,yp,turb_k_e,move_params,dt,space_cells)

f_phi = zeros(float_type, psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
write(base_filename*"array_shape", [psi_partions_num, psi_partions_num, y_res, x_res, nt])
if PSP_on 
    pptr.PSP_model!(f_phi ,xp, yp, turb_k_e, bc_interact, dt, "Uniform phi_1", mix_params, psi_mesh, space_cells, bc_params, verbose)
else
    pptr.make_f_phi_no_PSP!(f_phi,xp,yp, turb_k_e, bc_interact, "Uniform phi_1", psi_mesh, space_cells, bc_params, verbose)
end
write(base_filename*"data", f_phi)
println("data saved at:"*base_filename*"data")