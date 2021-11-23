import Random
include("PSP_Particletracking_module.jl")
pptr = PSPParticleTrackingReactions

Random.seed!(12345) #setting a seed

base_filename="Data/PSP_on_uniform_1_square_c_0_21_k_09_w_1_new_abs"

n=10
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
turb_k_e=0.938
bc_k = 0.25
PSP_on = true

space_cells = pptr.cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = pptr.psi_grid(psi_partions_num, psi_domain)

mix_params, move_params, bc_params = pptr.PSP_motion_bc_params(omega_bar, omega_sigma_2, C_0, B_format, c_phi, c_t, u_mean, bc_k, 10, true) #bc_params get dummy values

xp = zeros(np,nt+1)
yp = zeros(np,nt+1)
xp[:,1] = length_domain.*rand(np)
yp[:,1] = height_domain.*rand(np)
bc_interact = pptr.particle_motion_model(xp,yp,turb_k_e,move_params,dt,space_cells)

for NVP = [1:20;30;40;44:56;10 .*2 .^(3:10)]
    filename = base_filename*"_"*string(NVP)*"_vp"
    for bc_CLT = [true ]
        bc_CLT && (filename *= "_CLT")
        # isfile(filename) && (println(filename, ' ', "skipped");continue)

        bc_params = pptr.BC_params(bc_k, C_0, B_format, NVP, bc_CLT)
        f_phi = zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
        if PSP_on
            pptr.PSP_model!(f_phi,xp,yp, turb_k_e, bc_interact, dt, "Uniform phi_1", mix_params, psi_mesh, space_cells, bc_params, true)
        else
            pptr.make_f_phi_no_PSP!(f_phi,xp,yp, turb_k_e, bc_interact, "Uniform phi_1", psi_mesh, space_cells, bc_params, true)
        end
        write(filename, f_phi)
        println(filename, ' ', "success")
    end
end