import Random
include("PSP_Particletracking_module.jl")
pptr = PSPParticleTrackingReactions

Random.seed!(12345) #setting a seed

folders = ["Data","paper","PSP_on_uniform_1_5x100_02_001_c_0_21_k_1938_w_varied_new_abs_psi50_1000resparts_bck_03_T_02_ntmin","F32"]

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

T=0.2
length_domain = float_type(0.2) #length of periodic element
height_domain = float_type(0.01)
n=5
x_res=20*n #number of cells in x dim
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
PSP_on = true

nt = ceil(Int, T*2*sqrt((2/3)*turb_k_e)*min(x_res/length_domain,y_res/height_domain))   # number of time steps
dt = float_type(T/nt)#float_type(0.0002)  # time step
println(nt)
step_95=2*sqrt((2/3)*turb_k_e)*dt#95% of steps will be less than this
(step_95>length_domain/x_res || step_95>height_domain/y_res) && @warn "cell size smaller than >5% steps"

space_cells = pptr.cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = pptr.psi_grid(psi_partions_num, psi_domain)

_, move_params, bc_params_ = pptr.PSP_motion_bc_params(omega_bar, omega_sigma_2, C_0, B_format, c_phi, c_t, u_mean, bc_k, 10, true) #bc_params get dummy values

(bc_params_.bc_k*sqrt(pi*bc_params_.B/(bc_params_.C_0*turb_k_e))>1) && throw(ErrorException("reaction prob >1"))

xp = zeros(float_type, np,nt+1)
yp = zeros(float_type, np,nt+1)
xp[:,1] = length_domain.*rand(float_type, np)
yp[:,1] = height_domain.*rand(float_type, np)
bc_interact = pptr.particle_motion_model(xp,yp,turb_k_e,move_params,dt,space_cells)

for omega_bar = Float32.([1,4,8,12,16,18,20,22,24,28,32,40,50])
    mix_params = pptr.PSP_params(omega_bar, omega_sigma_2, c_phi, c_t)
    for NVP = [Float32(Inf)]
        (NVP < Inf) && (NVP=Int(NVP))
        filename = base_filename*string(NVP)*"_vp"
        for bc_CLT = [true ]
            bc_CLT && (filename *= "_CLT")
            filename *="_w"*string(omega_bar)
            PSP_on || (filename *= "_no_PSP")
            isfile(filename) && (println(filename, ' ', "skipped");continue)
            bc_params = pptr.BC_params(bc_k, C_0, B_format, NVP, bc_CLT)
            f_phi = zeros(float_type, psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
            if PSP_on
                pptr.PSP_model!(f_phi,xp,yp, turb_k_e, bc_interact, dt, "Uniform phi_1", mix_params, psi_mesh, space_cells, bc_params, true)
            else
                pptr.make_f_phi_no_PSP!(f_phi,xp,yp, turb_k_e, bc_interact, "Uniform phi_1", psi_mesh, space_cells, bc_params, true)
            end
            write(filename, f_phi)
            println(filename, ' ', "success")
        end
    end
end