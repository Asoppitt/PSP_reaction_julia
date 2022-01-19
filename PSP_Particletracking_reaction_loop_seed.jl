import Random
include("PSP_Particletracking_module.jl")
pptr = PSPParticleTrackingReactions


base_filename="Data/PSP_on_uniform_1_vlong_c_0_21_k_09_w_10_new_abs_Inf_vp_psi100_record_flux_u006"

n=5
x_res=20*n #number of cells in x dim
y_res=n
np = x_res*y_res*5000 # number of particles
dt = 0.002  # time step
nt = 60   # number of time steps
length_domain = 0.20 #length of periodic element
height_domain = 0.01
psi_partions_num = 100
psi_domain = [0.0,1.2]
omega_bar = 10.0
omega_sigma_2 = 0.25
C_0 = 1.0
B_format = "Constant"
c_phi = 25.0
c_t = 2.0
u_mean = 0.006
turb_k_e=0.938
bc_k = 0.25
PSP_on = true

space_cells = pptr.cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = pptr.psi_grid(psi_partions_num, psi_domain)

mix_params, move_params, bc_params = pptr.PSP_motion_bc_params(omega_bar, omega_sigma_2, C_0, B_format, c_phi, c_t, u_mean, bc_k, Inf, true)

xp = zeros(np,nt+1)
yp = zeros(np,nt+1)
f_phi = zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
flux = zeros(psi_partions_num, psi_partions_num, 1, x_res, nt)

N_runs = 5
seed_array = zeros(UInt32,4,10)
for i=1:N_runs
    rng=Random.seed!() #setting a seed
    seed_array[:,i] = rng.seed
    xp[:,1] = length_domain.*rand(np)
    yp[:,1] = height_domain.*rand(np)
    bc_interact = pptr.particle_motion_model_ref_start(xp,yp,turb_k_e,move_params,dt,space_cells)

    filename = base_filename*"/run"*string(i)
    if PSP_on
        pptr.PSP_model_inflow_record_flux!(f_phi, flux,xp,yp, turb_k_e, bc_interact, dt, "Uniform phi_1", mix_params, psi_mesh, space_cells, bc_params, true)
    else
        pptr.make_f_phi_no_PSP_inflow_record_flux!(f_phi, flux, xp,yp, turb_k_e, bc_interact, "Uniform phi_1", psi_mesh, space_cells, bc_params, true)
    end
    write(filename, flux)
    write(filename*"f_phi", f_phi)
    println(i, ' ', "success")
end
write(base_filename*"/_seed_array", seed_array)