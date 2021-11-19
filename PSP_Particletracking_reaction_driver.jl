import Random
include("PSP_Particltracking_module.jl")
pptr = PSPParticleTrackingReactions

Random.seed!(12345) #setting a seed

filename="Data/PSP_on_uniform_1_smolsquare_c_0_21_k_09_w_1_new_abs_80_vp_CLT"

n=5
x_res=n #number of cells in x dim
y_res=n
np = x_res*y_res*1000 # number of particles
dt = 0.005  # time step
nt = 60   # number of time steps
length_domain = 0.10 #length of periodic element
height_domain = 0.04
psi_partions_num = 20
psi_domain = [0.0,1.2]
omega_bar = 1.0
omega_sigma_2 = 0.25
C_0 = 2.1
B_format = "Constant"
c_phi = 25.0
c_t = 2.0
u_mean = 0.0
turb_e_init=0.938
bc_k = 0.25
num_vp = 80
bc_CLT = true
PSP_on = true

space_cells = pptr.cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = pptr.psi_grid(psi_partions_num, psi_domain)
mix_params, move_params, bc_params = pptr.PSP_motion_bc_params(omega_bar, omega_sigma_2, C_0, B_format, c_phi, c_t, u_mean, bc_k, num_vp, bc_CLT)

xp = zeros(np,nt+1)
yp = zeros(np,nt+1)
xp[:,1] = length_domain.*rand(np)
yp[:,1] = height_domain.*rand(np)
bc_interact = pptr.particle_motion_model(xp,yp,turb_k_e,move_params,dt,space_cells)

f_phi = zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
if PSP_on
    pptr.PSP_model!(f_phi,xp,yp, turb_k_e, bc_interact, dt, "Uniform phi_1", mix_params, psi_mesh, space_cells, bc_params)
else
    pptr.make_f_phi_no_PSP!(f_phi,xp,yp, turb_k_e, bc_interact, "Uniform phi_1", psi_mesh, space_cells, bc_params)
end
write(filename, f_phi)