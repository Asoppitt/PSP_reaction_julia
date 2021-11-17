using Random, Distributions
import StatsBase as sb
import LinearAlgebra as la
import Statistics as st

Random.seed!(12345) #setting a seed
Initial_condition = "Uniform phi_1"
phi_eps = 10^(-8) #small non-zero phi, used to replace 0, as histogram requires phi>0
psi_partions_num = 20 #number of partitions in 1 state space direction
phi_domain = [0,1.2]
psi_1 = LinRange(phi_domain[1], phi_domain[2], psi_partions_num+1) #defines boundaries for cells in psi space
psi_2 = LinRange(phi_domain[1], phi_domain[2], psi_partions_num+1)
CLT_limt = true

function assign_pm!(phi_pm::Matrix{Int}, phi_array_t::Array{Float64}, particles::Vector{Int}, cell_points::Vector{Int}) 
    # might be worth a strategy that removes tested pairs, can't see how to make it not require a biiiig temp variable though 
    for particle in particles
        while true
            try_pm = sb.sample(cell_points,2,replace=false, ordered=true) #ordering added to reduce sample space as order for pm pairs is irrlevent
            test_dot = la.dot((phi_array_t[:,try_pm[1]]-phi_array_t[:,particle]),(phi_array_t[:,try_pm[2]]-phi_array_t[:,particle]))
            if  test_dot<= 0
                phi_pm[:,particle] = try_pm
                break
            end
        end
    end
end

function assign_pm_single!(phi_pm::Matrix{Int}, phi_array_t::Array{Float64}, particles::Vector{Int}, cell_points::Vector{Int}, index_tbc::Int)
    known_index = 1-(index_tbc-1)+1 #flips 1 to 2 and 2 to 1
    for particle in particles
        !in(particle, cell_points) && throw(ArgumentError("Particles must be in cell"))
        cell_points_set = Set(cell_points)
        while true
            try_pm = rand(cell_points_set) 
            test_dot = la.dot((phi_array_t[:,try_pm[1]]-phi_array_t[:,particle]),(phi_array_t[:,phi_pm[known_index,particle]]-phi_array_t[:,particle]))
            if  test_dot<= 0
                phi_pm[index_tbc,particle] = try_pm[1]
                break
            end
            setdiff!(cell_points_set,try_pm)
        end
    end
end
if Initial_condition == "Uniform phi_1"
    function set_phi_as_ic!(phi_array::Array{Float64},yp::Array{Float64})
        nparticles = size(phi_array)[2]
        phi_array[2,:] = abs.(phi_eps*randn(nparticles)) #pdf can't find zeros
        phi_array[1,:] .= 1;
    end
elseif Initial_condition == "triple delta"
    function set_phi_as_ic!(phi_array::Array{Float64},yp::Array{Float64})
        nparticles = size(phi_array)[2]
        local delta_selector = rand(1:3, nparticles)
        local noise_term = randn(nparticles)

        phi_array[1,delta_selector.=1] = -sqrt(3)/2 .+phi_eps .*noise_term[delta_selector.=1]
        phi_array[2,delta_selector.=1] .= -0.5

        phi_array[1,delta_selector.=2] = sqrt(3)/2 .+phi_eps .*noise_term[delta_selector.=2]
        phi_array[2,delta_selector.=2] .= -0.5

        phi_array[1,delta_selector.=3] .= phi_eps
        phi_array[2,delta_selector.=3] = 1.0 .+phi_eps .*noise_term[delta_selector.=3]
    end
elseif Initial_condition == "2 layers"
    function set_phi_as_ic!(phi_array::Array{Float64},yp::Array{Float64})
        nparticles = size(phi_array)[2]
        local noise_term = randn(nparticles)
        # local uniform_noise = rand(nparticles).-0.5

        phi_array[2,[yp.>0.5*height_domain]] = abs.(phi_eps*noise_term[yp.>0.5*height_domain] )
        phi_array[1,[yp.>0.5*height_domain]] .= 1

        phi_array[1,[yp.<=0.5*height_domain]] = abs.(phi_eps*noise_term[yp.<=0.5*height_domain] )
        phi_array[2,[yp.<=0.5*height_domain]] .= 1 #.+ uniform_noise[yp[particles,1].<=0.5*height_domain].*0.05
    end
elseif Initial_condition == "double delta"
    function set_phi_as_ic!(phi_array::Array{Float64},yp::Array{Float64})
        nparticles = size(phi_array)[2]
        local noise_term = randn(nparticles)
        local delta_selector = rand([true,false], nparticles)
        # local uniform_noise = rand(nparticles).-0.5

        phi_array[2,delta_selector] = abs.(phi_eps*noise_term[delta_selector] )
        phi_array[1,delta_selector] .= 1 #.+ uniform_noise[delta_selector.==1].*0.05

        phi_array[1,delta_selector] = abs.(phi_eps*noise_term[.!delta_selector] )
        phi_array[2,delta_selector] .= 1 #.+ uniform_noise[delta_selector.==2].*0.05
    end
end

function assign_f_phi_cell!(f_phi_cell::Array{Float64},phi_array::Array{Float64})
    "use if cell_points has alreday been determined"
    for psi1_i=1:psi_partions_num
        in_1 = psi_1[psi1_i].<=phi_array[1,:].<psi_1[psi1_i+1]
        for psi2_i=1:psi_partions_num
            in_2 = psi_2[psi2_i].<=phi_array[2,:].<psi_2[psi2_i+1]
            f_phi_cell[:,:] .= sum(in_1.&in_2)
        end
    end
    f_phi_cell[:, :] = f_phi_cell./(size(phi_array)[2])
end 

function assign_f_phi!(f_phi_t::Array{Float64},phi_array::Array{Float64})
    for i in 1:y_res
        in_y = y_edges[i].<=yp[:,t_index].<y_edges[i+1]
        for j in 1:x_res
            in_x = x_edges[j].<=xp[:,t_index].<x_edges[j+1]
            cell_particles = findall(in_x.&in_y)
            assign_f_phi_cell!(f_phi_t[:,:,i, j],phi_array[:,cell_particles])
        end 
    end
end 
if CLT_limt
    function bc_absorbtion(abs_points::Array{TF}, nvpart_per_part::Int, turb_k_e::Array{TF}) where TF<:AbstractFloat
        n_abs = size(abs_points)[2]
        abs_k = bc_k.*ones(2,n_abs)
        effective_v_particles =( abs_points.*nvpart_per_part)
        #K for Erban and Chapman approximation 
        P = zeros(2,n_abs)
        P[1,:] = min.(abs_k[1,:].*sqrt.(B.*pi./(C_0.*turb_k_e)),1)
        P[2,:] = min.(abs_k[2,:].*sqrt.(B.*pi./(C_0.*turb_k_e)),1)
        #by CLT approx dist for number of virtual particles to have reacted
        xi = randn(2,n_abs).*sqrt.((P.*(1 .-P)))
        #catching places where all mass has been removed
        xi = [effective_v_particles[i,j]>0 ? xi[i,j]./sqrt(effective_v_particles[i,j]) : 0 for i in 1:2, j in 1:n_abs]
        ratios = max.(min.((1 .-P) + xi,1),0)
        abs_points[:, :] = abs_points.*ratios
    end
else
    #using binomal noise for small numbers of vparticles
    function bc_absorbtion(abs_points::Array{TF}, nvpart_per_part::Int, turb_k_e::Array{TF}) where TF<:AbstractFloat
        n_abs = size(abs_points)[2]
        abs_k = bc_k.*ones(2,n_abs)
        effective_v_particles =( abs_points.*nvpart_per_part)
        #K for Erban and Chapman approximation 
        P = zeros(2,n_abs)
        P[1,:] = min.(abs_k[1,:].*sqrt.(B.*pi./(C_0.*turb_k_e)),1)
        P[2,:] = min.(abs_k[2,:].*sqrt.(B.*pi./(C_0.*turb_k_e)),1)
        #Binomal dist for number of virtual particles to have reacted
        xi_dist = Binomial.(ceil.(effective_v_particles),1 .-P)
        xi = [rand(xi_dist[i,j]) for i in 1:2, j in 1:n_abs]
        ratios = [effective_v_particles[i,j]>0 ? xi[i,j]./ceil.(effective_v_particles[i,j]) : 0 for i in 1:2, j in 1:n_abs]
        abs_points[:, :] = abs_points.*ratios
    end
end
function bc_absorbtion(abs_points::Array{TF}, nvpart_per_part::AbstractFloat, turb_k_e::Array{TF}) where TF<:AbstractFloat
    nvpart_per_part != Inf && throw(DomainError("nvpart_per_part must be Int or Inf"))
    n_abs = size(abs_points)[2]
    abs_k = bc_k.*ones(2,n_abs)
    #K for Erban and Chapman approximation 
    P = zeros(2,n_abs)
    P[1,:] = min.(abs_k[1,:].*sqrt.(B.*pi./(C_0.*turb_k_e)),1)
    P[2,:] = min.(abs_k[2,:].*sqrt.(B.*pi./(C_0.*turb_k_e)),1)
    ratios = 1 .-P #taking mean for limiting case
    abs_points[:, :] = abs_points.*ratios
end
for test_index=1:1
filename="Data/PSP_on_uniform_1_smolsquare_c_0_21_k_09_w_1_new_abs_80_vp_CLT"

n=5 #a convinent method of increascing resolution while maintaining 
    # resolution ratios and number of particles per cell

x_res=n #number of cells in x dim
y_res=n
np = x_res*y_res*1000 # number of particles
dt = 0.005  # time step
nt = 60   # number of time steps
length_domain = 0.10 #length of periodic element
height_domain = 0.04
omega_shape_init = 0.25
omega_sigma_2 = 0.25
c_phi = 25
c_t = 2
u_max = 0
function u_mean(y)
    u_max.*ones(np)#-4*u_max.*y.*(y.-height_domain)./height_domain^2
end
C_0 = 2.1  # rate parameter for velocity change
B=(1.5*C_0) #for reducing the langevin to a standard form - part of boundary conditions implementaion from Erban and Chapman 06
PSP_off = true
new_inflow = false
record_BC_flux = true
#assuming each psp particle consists of a number of virtual particles with binary scalar values 
#when bc is computed, the distribution absorbed is based on this number
nvpart_per_part = 80

bc_k=0.25
omega_mean=1
turb_e_init=0.938
T_omega = 1/(omega_mean); #approximation to frequency timescale
phip = zeros((2, np, nt+1)) #scalar concentration at these points
#ititialising position
xp = length_domain*rand(Float64, (np,nt+1)) # position of particle in x-direction
yp = height_domain*rand(Float64, (np,nt+1)) # position of paticle in y-direction


record_BC_flux && (flux_y0 = zeros(2,x_res,nt))
record_BC_flux && (flux_y0_vel = zeros(x_res,nt))
dphi_x_phi = ((phi_domain[2]-phi_domain[1])/psi_partions_num)^2 #approx for the phi-space 2-D lebesgue measure that f_phi is defined wrt

set_phi_as_ic!(phip[:,:,1],yp[:,1])

omegap = zeros(np,nt+1) #turbulence frequency
omega0_dist = Gamma(omega_mean/omega_shape_init,omega_shape_init)
omegap[:,1] = rand(omega0_dist, np)

x_edges = LinRange(0,length_domain,x_res+1)
y_edges = LinRange(0,height_domain,y_res+1)

#array for storing indicies of paired particles for psp
phi_pm = zeros(Int, 2, np)

f_phi = zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1) # histogram of each cell at a given time 



#initialising phi_pm and f_phi
for i in 1:y_res
    in_y = y_edges[i].<yp[:,1].<y_edges[i+1]
    for j in 1:x_res
        in_x = x_edges[j].<xp[:,1].<x_edges[j+1]
        cell_particles = findall(in_x.&in_y)
        assign_pm!(phi_pm, phip[:,:,1], cell_particles, cell_particles)
        assign_f_phi_cell!(f_phi[:,:,i,j,1],phip[:,cell_particles,1])
    end 
end

#time stamp until new p/m bound found, needed to ensure particles are
#decorreltaed
t_decorr_p = 1 ./(c_t.*omegap[phi_pm[1,:],1])
t_decorr_m = 1 ./(c_t.*omegap[phi_pm[2,:],1])

k= ones(np).*turb_e_init

#intitial vaules of velocity, maintainingconsitancy with energy
uxp = randn(np).*sqrt.(2/3*k)
uyp = randn(np).*sqrt.(2/3*k)

for t in 1:nt
    println(t)
    #Using the SLM model for particle location
    for i in 1:y_res
        in_y = y_edges[i].<yp[:,t].<y_edges[i+1]
        for j in 1:x_res
            in_x = x_edges[j].<xp[:,t].<x_edges[j+1]
            cell_particles = findall(in_x.&in_y)
            k[cell_particles].=turb_e_init#0.5*(st.mean(uxp[cell_particles].^2)+st.mean(uyp[cell_particles].^2))*1.5 #turb_e_init;
        end
    end
    uxp = uxp+(-0.5*B*omega_mean*uxp)*dt+randn(np).*sqrt.(C_0.*k.*omega_mean.*dt); 
    uyp = uyp+(-0.5*B*omega_mean*uyp)*dt+randn(np).*sqrt.(C_0.*k.*omega_mean.*dt); 
    uxp_full = u_mean(yp[:,t])+uxp
    xp[:,t+1] = xp[:,t] + (uxp_full)*dt # random walk in x-direction
    yp[:,t+1] = yp[:,t] + uyp*dt # random walk in y-direction


    # Reflection particles at boundaries

    # Reflection at upper boundary y>height_domain
    # doing closed on top open on bottom, as cell detection is open on top,
    # closed on bottom
    mag = findall(yp[:,t+1].>=height_domain) # index of particle with yp>height_domain
    dim_mag = size(mag) # dimension of array "mag"

    y_mag_succ = yp[mag,t+1] # yp at time t+1 corresponding to the index "mag"

    V1 = height_domain.*ones(dim_mag) 

    ypr_mag = V1*2 .- y_mag_succ  # yp at time t+1 of the reflected particle

    yp[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
    uyp[mag] = -uyp[mag] #reflecting velocity

    # bc_absorbtion(phip[:,mag,t+1], nvpart_per_part, k[mag]) disabled to match paper

    # Reflection at lower boundary y<0
    mag = findall(yp[:,t+1].<=0) # index of particle with yp>height_domain
    dim_mag = size(mag) # dimension of array "mag"

    record_BC_flux && (missing_mass = phip[:,:,t])
    bc_absorbtion(phip[:,mag,t+1], nvpart_per_part, k[mag])
    if record_BC_flux
        missing_mass -= phip[:,:,t]
        in_y = y_edges[1].<yp[:,t].<y_edges[2] #is in lowest row?
        for j in 1:x_res
            in_x = x_edges[j].<xp[:,t].<x_edges[j+1]
            cell_particles = findall(in_x.&in_y)
            flux_y0[:,j,t] = mean(missing_mass[:,cell_particles],dims=2)[:,1]
            flux_y0_vel[j,t] = mean(uyp[mag])
        end
    end
    y_mag_succ = yp[mag,t+1] # yp at time t+1 corresponding to the index "mag"

    ypr_mag = - y_mag_succ  # yp at time t+1 of the reflected particle

    yp[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
    uyp[mag] = -uyp[mag] #reflecting velocity


    #bc at end (y=length_domain) of domain
    end_indicies = xp[:,t+1].>=length_domain #index of particle with xp>length
    dim_end = sum(end_indicies)

    end_x = xp[end_indicies,t+1]
    xpr_end = end_x .- length_domain #shifting particles back to begining
    xp[end_indicies,t+1] = xpr_end #replacing x coords

    if new_inflow
        #resetting concentration of "new" particles
        set_phi_as_ic!(phip[:,end_indicies,1],yp[end_indicies,1])

        #resetting energy of "new" particles
        uxp[end_indicies] = randn(dim_end,1)*sqrt(2/3*turb_e_init)
        uyp[end_indicies] = randn(dim_end,1)*sqrt(2/3*turb_e_init)
        k[end_indicies] .= turb_e_init
    end

    #bc at start (y=0) of domain
    start_indicies = xp[:,t+1].<=0 #index of particle with xp>length
    dim_start = sum(start_indicies)

    if new_inflow
        xpr_start = -xp[start_indicies,t+1] #reflection
        xp[start_indicies,t+1] = xpr_start #replacing x coords
        uxp[start_indicies] = - uxp[start_indicies]

        #resetting concentration of "new" particles
        set_phi_as_ic!(phip[:,start_indicies,1],yp[start_indicies,1])

        #resetting energy of "new" particles
        uxp[start_indicies] = randn(dim_start,1)*sqrt(2/3*turb_e_init)
        uyp[start_indicies] = randn(dim_start,1)*sqrt(2/3*turb_e_init)
        k[start_indicies] .= turb_e_init
    else
        #true periodic bc
        xpr_start = length_domain .+ xp[start_indicies,t+1] 
        xp[start_indicies,t+1] = xpr_start #replacing x coords
    end


    #PSP begins here
    if PSP_off
        dphi = zeros(2,np)
    else
        #E-M solver for omega 
        dw = sqrt(dt).*randn(np) #random draws
        omegap[:,t+1] = omegap[:,t]-(omegap[:,t].-omega_mean)./T_omega.*dt + sqrt.(omegap[:,t].*(2*omega_sigma_2*omega_mean/T_omega)).*dw
        omegap[:,t+1] = omegap[:,t+1].*(omegap[:,t+1].>0) #enforcing positivity

        #stepping the decorrelation times
        t_decorr_p = t_decorr_p.-dt;
        t_decorr_m = t_decorr_m.-dt;
        #sets index of those to be renewed with 0 - which doesn't correspond to any particle
        phi_pm[1,(t_decorr_p.<=0)] .= 0 #the t_decorrs are 2d for some reason
        phi_pm[2,(t_decorr_m.<=0)] .= 0

        #split into cells, compute centres/targets, run ODE step
        for i in 1:y_res
            in_y = y_edges[i].<=yp[:,t+1].<y_edges[i+1]
            for j in 1:x_res
                in_x = x_edges[j].<=xp[:,t+1].<x_edges[j+1]
                cell_particles = findall(in_x.&in_y)
                if 0 in cell_particles
                    print("0 in cell_particles ", i,j,' ',t,'\n')
                end
                (length(cell_particles)==0) && throw(BoundsError(cell_particles))
                #update pairs to ensure all are within the same bounds
                p_nin = .!(in.(phi_pm[1,cell_particles],Ref(cell_particles)))
                m_nin = .!(in.(phi_pm[2,cell_particles],Ref(cell_particles)))
                pm_nin = p_nin .& m_nin
                p_nin = xor.(p_nin , pm_nin)
                m_nin = xor.(m_nin , pm_nin)
                (sum(p_nin)>0)&&assign_pm_single!(phi_pm, phip[:,:,t],cell_particles[p_nin], cell_particles, 1)
                (sum(m_nin)>0)&&assign_pm_single!(phi_pm, phip[:,:,t],cell_particles[m_nin], cell_particles, 2)
                (sum(pm_nin)>0)&&assign_pm!(phi_pm, phip[:,:,t], cell_particles[pm_nin], cell_particles)
                t_decorr_p[cell_particles[pm_nin.&p_nin]] = 1 ./(c_t.*omegap[phi_pm[1,cell_particles[pm_nin.&p_nin]],1])
                t_decorr_m[cell_particles[pm_nin.&p_nin]] = 1 ./(c_t.*omegap[phi_pm[2,cell_particles[pm_nin.&m_nin]],2])
            end 
        end
        phi_c = 0.5.*(phip[:,phi_pm[1,:],t]+phip[:,phi_pm[2,:],t])
        diffusion = zeros(2,np)
        diffusion[1,:] = c_phi.*0.5.*omegap[:,t].*(phip[1,:,t]-phi_c[1,:])
        diffusion[2,:] = c_phi.*0.5.*omegap[:,t].*(phip[2,:,t]-phi_c[2,:])
        reaction = zeros(2,np) # body reaction
        dphi = (-diffusion .+ reaction).*dt

        # ensuring mean 0 change
        # generating a random orthonormal basis
        # is 2-d so genrate a random unit vector from an angle and proceed based
        # on that
        angle = 2*pi*rand(1)[1];
        e_1 = [cos(angle),sin(angle)]
        handedness = sb.sample([-1,1],1)[1] #randomly choose betwen left or right handed system
        e_2 = handedness*[e_1[2],-e_1[1]]
        T=zeros(2,2)
        T[:,1] = e_1  #coord transform matrix
        T[:,2] = e_2
        dphi = T\dphi  #transform to new coords
        #performing adjustment to mean 0
        corr_factor = zeros(2,np)
        for i in 1:y_res
            in_y = y_edges[i].<=yp[:,t+1].<y_edges[i+1]
            for j in 1:x_res
                in_x = x_edges[j].<=xp[:,t+1].<x_edges[j+1]
                cell_particles = findall(in_x.&in_y)
                for phi_i=1:2
                    phi_mean = mean(dphi[phi_i,cell_particles])
                    if phi_mean != 0
                        cell_points_pos = dphi[phi_i,cell_particles].>0
                        cell_points_neg = dphi[phi_i,cell_particles].<0
                        phi_pos_mean = mean(cell_points_pos.*dphi[phi_i,cell_particles])
                        phi_neg_mean = mean((cell_points_neg).*dphi[phi_i,cell_particles])
                        if phi_mean>0
                            corr_factor[phi_i,cell_particles]=.-cell_points_pos*(phi_neg_mean./phi_pos_mean) + (1 .- cell_points_pos)
                            (any(isinf.(corr_factor[phi_i,cell_particles])) || any(isnan.(corr_factor[phi_i,cell_particles]))) && print(" +corr ")
                        else
                            corr_factor[phi_i,cell_particles]=.-(cell_points_neg)*(phi_pos_mean./phi_neg_mean) + (1 .- cell_points_neg)
                            (any(isinf.(corr_factor[phi_i,cell_particles])) || any(isnan.(corr_factor[phi_i,cell_particles]))) && print(" -corr ")
                        end
                    end
                end
            end
        end
        dphi = corr_factor.*dphi
        dphi = T*dphi #return to old coords
    end

    phip[:,:,t+1] = phip[:,:,t]+dphi
    if !(Initial_condition == "triple delta")
        phip[:,:,t+1] = phip[:,:,t+1].*(phip[:,:,t+1].>0) #forcing positive concentration
    end

    assign_f_phi(f_phi[:,:,:,:,t],phip[:,:,t])

end

write(filename, f_phi)
record_BC_flux && write(filename*"flux", flux_y0)
record_BC_flux && write(filename*"flux_vel", flux_y0_vel)
print("Success")
end