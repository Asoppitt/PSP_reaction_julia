import Random
using TurbulentMixingParticleTrackingReactions


base_filename="Data/paper/PSP_off_uniform_1_10x50_005_001_c_0_21_k_0938_w_varied_new_abs_Inf_vp_psi50_2000resparts_bck_03_dt0001"

# run(`mkdir $base_filename`)

n=10
x_res=n*5 #number of cells in x dim
y_res=n
np = x_res*y_res*2000 # number of particles
dt = 0.001  # time step
nt = 200   # number of time steps
length_domain = 0.05 #length of periodic element
height_domain = 0.01
psi_partions_num = 50
psi_domain = [0.0,1.2]
omega_bar = [22.5,25]
# append!(omega_bar,2.5:2.5:25)
omega_sigma_2 = 0.25
C_0 = 2.1
B_format = "Constant"
c_phi = 25.0
c_t = 2.0
u_mean = 0.006
turb_k_e=0.938
bc_k=0.3#(x) = x/(x+0.5)
PSP_on = false

space_cells = cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = psi_grid(psi_partions_num, psi_domain)


xp = zeros(np,nt+1)
yp = zeros(np,nt+1)
f_phi = zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
flux = zeros(psi_partions_num, psi_partions_num, 1, x_res, nt)

write(base_filename*"/test",' ')#check folders exist

N_runs = length(omega_bar)
seed_array = zeros(UInt32,4,N_runs)
for i=1:N_runs
    mix_params, move_params, bc_params = PSP_motion_bc_params(omega_bar[i], omega_sigma_2, C_0, B_format, c_phi, c_t, u_mean, bc_k, Inf, true)
    rng=Random.seed!() #setting a seed
    seed_array[:,i] = rng.seed
    xp[:,1] = length_domain.*rand(np)
    yp[:,1] = height_domain.*rand(np)
    bc_interact = particle_motion_model_ref_start(xp,yp,turb_k_e,move_params,dt,space_cells)

    filename = base_filename*"/omega"*string(omega_bar[i])
    if PSP_on
        PSP_model!(f_phi, xp,yp, turb_k_e, bc_interact, dt, "Uniform phi_1", mix_params, psi_mesh, space_cells, bc_params, true)
    else
        make_f_phi_no_PSP!(f_phi, xp,yp, turb_k_e, bc_interact, "Uniform phi_1", psi_mesh, space_cells, bc_params, true)
    end
    write(filename, f_phi)
    println(i, ' ', "success")
end
write(base_filename*"/_seed_array", seed_array)