import Random
using TurbulentMixingParticleTrackingReactions
# include("/home/u2093579/Documents/Project_code/TurbulentMixingParticleTrackingReactions/src/PSP_Particletracking_module.jl")

Random.seed!(12345) #setting a seed

folders = ["Data","test_DNS_Params","F32"]

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
length_domain = float_type(5e-4) #length of periodic element
height_domain = float_type(2/5*5e-4)
n=5
x_res=5*n #number of cells in x dim
y_res=2*n
np = x_res*y_res*1000 # number of particles
psi_partions_num = 50
psi_domain = [float_type(0.0),float_type(1.01)]
omega_bar = float_type(3642.765106925031)
omega_sigma_2 = float_type(omega_bar*1000)
C_0 = float_type(2.1)
B_format = "Constant"
c_phi = float_type(25.0)
c_t = float_type(2.0)
u_mean = float_type(0.01)
turb_k_e= float_type(0.00020448453573087953)#float_type(0.938)
bc_k = float_type(0.6)#for langmuir,is max sorbtion 
# K=float_type(0.25)#0.125
lang_half_sat=float_type(0.125)#langmuir saturation constant, value of phi for which removal rate is halfed
bc_non_lin_corr(phi) = lang_half_sat ./(lang_half_sat .+phi)
# bc_non_lin_corr=((c1,c2)->-5*(c2-K/c1), (c1,c2)->-5*(c1-K/c2))
PSP_on = true
ic = ("Uniform phi_1")
# reaction = ((c1,c2)->-10*(c1*c2-K), (c1,c2)->-10*(c1*c2-K))
reaction = ((c1,c2)->0.0, (c1,c2)->0.0)

u_99_f =3*sqrt((2/3)*turb_k_e)
u_95_f =2*sqrt((2/3)*turb_k_e)
nt = max(100,ceil(Int, T*max((x_res/length_domain)*(u_99_f+u_mean),(y_res/height_domain)*2*u_99_f)))   # number of time steps
dt = float_type(T/nt)#float_type(0.0002)  # time step
println(nt)
#95% of steps will be less than this
(((u_95_f+u_mean)*dt)>length_domain/x_res || ((u_95_f)*dt)>height_domain/y_res) && @warn "cell size smaller than >5% steps"
space_cells = cell_grid(x_res,y_res,length_domain,height_domain)
psi_mesh = psi_grid(psi_partions_num, psi_domain)

# _, move_params, bc_params_ = PSP_motion_bc_params(omega_bar, omega_sigma_2, C_0, B_format, c_phi, c_t, u_mean, bc_k, 10, true,bc_non_lin_corr) #bc_params get dummy values

# (bc_params_.bc_k*sqrt(pi*bc_params_.B/(bc_params_.C_0*turb_k_e))>1) && throw(ErrorException("reaction prob >1"))

# xp = zeros(float_type, np,nt+1)
# yp = zeros(float_type, np,nt+1)
# xp[:,1] = length_domain.*rand(float_type, np)
# yp[:,1] = height_domain.*rand(float_type, np)
# bc_interact = particle_motion_model(xp,yp,turb_k_e,move_params,dt,space_cells)
# any(xp.>length_domain)||any(xp.<0) && throw(ErrorException(""))
# throw(ErrorException("n"))
write(base_filename*"array_shape", [psi_partions_num, psi_partions_num, y_res, x_res, nt])
for omega_bar = float_type.([omega_bar])
    # omega_sigma_2 = float_type(0.25*omega_bar)
    
    mix_params, move_params, _ =PSP_motion_bc_params(omega_bar, omega_sigma_2,float_type(0.0), C_0, B_format, c_phi, c_t, u_mean, bc_k, 200, corr_func=bc_non_lin_corr, bulk_reaction=reaction)
    xp = zeros(float_type, np,nt+1)
    yp = zeros(float_type, np,nt+1)
    xp[:,1] = length_domain.*rand(float_type, np)
    yp[:,1] = height_domain.*rand(float_type, np)
    bc_interact = particle_motion_model(xp,yp,turb_k_e,move_params,dt,space_cells)
    any(xp.>length_domain)||any(xp.<0) && throw(ErrorException(""))
    for NVP = [float_type(Inf)]
        (NVP < Inf) && (NVP=Int(NVP))
        filename = base_filename*string(NVP)*"_vp"
        for bc_CLT = [true ]
            bc_CLT && (filename *= "_CLT")
            filename *="_w"*string(omega_bar)
            PSP_on || (filename *= "_no_PSP")
            # isfile(filename) && (println(filename, ' ', "skipped");continue)
            bc_params = BC_params(bc_k, C_0, B_format, NVP, corr_func=bc_non_lin_corr)
            f_phi = zeros(float_type, psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
            # phigamma_mean = zeros(float_type, y_res, x_res, nt)
            edge_mean = zeros(float_type,nt+1)
            edge_2 = zeros(float_type,nt+1)
            edge_2_v = zeros(float_type,nt)
            if PSP_on
                # PSP_model!(f_phi,xp,yp, turb_k_e, bc_interact, dt, ic, mix_params, psi_mesh, space_cells, bc_params, true)
                # PSP_model_record_phi_local_diff!(phigamma_mean ,f_phi ,xp,yp, turb_k_e, bc_interact, dt, ic, mix_params, psi_mesh, space_cells, bc_params, true)
                PSP_model_record_reacting_mass!(edge_mean, edge_2, edge_2_v, f_phi,xp, yp, turb_k_e, bc_interact, dt, ic,  mix_params, psi_mesh, space_cells, bc_params, true)
            else
                make_f_phi_no_PSP!(f_phi,xp,yp, turb_k_e, bc_interact, ic, psi_mesh, space_cells, bc_params, true)
            end
            write(filename, f_phi)
            write(filename*"edge_mean", edge_mean)
            write(filename*"edge_squared", edge_2)
            write(filename*"edge_squared_velocity", edge_2_v)
            # write(filename*"phigamma", phigamma_mean)
            println(filename, ' ', "success")
        end
    end
end