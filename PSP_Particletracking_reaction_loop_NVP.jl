import Random
include("PSP_Particletracking_module.jl")
pptr = PSPParticleTrackingReactions

Random.seed!(12345) #setting a seed

base_filename="Data/f32test/PSP_on_uniform_1_10x50_02_004_c_0_21_k_0938_w_varied_new_abs_Inf_vp_psi50_2000resparts_bck_025_dt0001"

run(`mkdir $base_filename`)

float_type = Float32
n=10
x_res=5*n #number of cells in x dim
y_res=n
np = x_res*y_res*2000 # number of particles
dt = float_type(0.001)  # time step
nt = 200   # number of time steps
length_domain = float_type(0.20) #length of periodic element
height_domain = float_type(0.04)
psi_partions_num = 50
psi_domain = [float_type(0.0),float_type(1.2)]
omega_bar = float_type(5.0)
omega_sigma_2 = float_type(0.25)
C_0 = float_type(2.1)
B_format = "Constant"
c_phi = float_type(25.0)
c_t = float_type(2.0)
u_mean = float_type(0.1)
turb_k_e= float_type(0.938)
bc_k = float_type(0.25)
PSP_on = true

space_cells = pptr.cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = pptr.psi_grid(psi_partions_num, psi_domain)

mix_params, move_params, bc_params = pptr.PSP_motion_bc_params(omega_bar, omega_sigma_2, C_0, B_format, c_phi, c_t, u_mean, bc_k, 10, true) #bc_params get dummy values

(bc_params.bc_k*sqrt(pi*bc_params.B/(bc_params.C_0*turb_k_e))>1) && throw(ErrorException("reaction prob >1"))

xp = zeros(float_type, np,nt+1)
yp = zeros(float_type, np,nt+1)
xp[:,1] = length_domain.*rand(float_type, np)
yp[:,1] = height_domain.*rand(float_type, np)
bc_interact = pptr.particle_motion_model(xp,yp,turb_k_e,move_params,dt,space_cells)

for NVP = [20,50,100,200,Inf]
    (NVP < Inf) && (NVP=Int(NVP))
    filename = base_filename*"/"*string(NVP)*"_vp"
    for bc_CLT = [true ]
        bc_CLT && (filename *= "_CLT")
        isfile(filename) && (println(filename, ' ', "skipped");continue)
        bc_params = pptr.BC_params(bc_k, C_0, B_format, NVP, bc_CLT)
        f_phi = zeros(float_type, psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
        if PSP_on
            @time pptr.PSP_model!(f_phi,xp,yp, turb_k_e, bc_interact, dt, "Uniform phi_1", mix_params, psi_mesh, space_cells, bc_params, true)
        else
            @time pptr.make_f_phi_no_PSP!(f_phi,xp,yp, turb_k_e, bc_interact, "Uniform phi_1", psi_mesh, space_cells, bc_params, true)
        end
        write(filename, f_phi)
        println(filename, ' ', "success")
    end
end