struct CellGrid{T<:AbstractFloat}
    x_res ::Int
    y_res ::Int
    length_domain ::T
    height_domain ::T
    x_edges ::AbstractVector{T}
    y_edges ::AbstractVector{T}
end

struct PsiGrid{T<:AbstractFloat}
    psi_partions_num ::Int
    phi_domain ::AbstractVector{T}
    psi_1 ::AbstractVector{T}
    psi_2 ::AbstractVector{T}
end

struct PSPParams{T<:AbstractFloat}
    omega_bar::T
    omega_sigma_2::T
    T_omega::T
    C_0::T
    B::T
    c_phi::T
    c_t::T
end

struct MotionParams{T<:AbstractFloat}
    omega_bar::T
    C_0::T
    B::T
    u_mean::T
end

function cell_grid(x_res ::Int,y_res ::Int,length_domain ::T,height_domain ::T) where T<:AbstractFloat #constructor for CellGrid
    return CellGrid(x_res,y_res,length_domain,height_domain,LinRange(0,length_domain,x_res+1),LinRange(0,height_domain,y_res+1))
end

function psi_grid(psi_partions_num ::Int, phi_domain ::AbstractVector{T}) where T<:AbstractFloat #constructor for PsiGrid
    return PsiGrid(psi_partions_num, phi_domain,LinRange(phi_domain[1], phi_domain[2], psi_partions_num+1),LinRange(phi_domain[1], phi_domain[2], psi_partions_num+1))
end

function PSP_params(omega_bar::T, omega_sigma_2::T, C_0::T, B_format::String, c_phi::T, c_t::T) where T<:AbstractFloat
    if B_format == "Decay"
        return PSPParams(omega_bar, omega_sigma_2, 1/omega_bar, C_0, (1+1.5*C_0), c_phi, c_t)
    elseif B_format == "Constant"
        return PSPParams(omega_bar, omega_sigma_2, 1/omega_bar, C_0, (1.5*C_0), c_phi, c_t)
end

function motion_params(omega_bar::T,C_0::T, B_format::String, u_mean::T) where T<:AbstractFloat
if B_format == "Decay"
    return MotionParams(omega_bar, C_0, (1+1.5*C_0), u_mean)
elseif B_format == "Constant"
    return MotionParams(omega_bar, C_0, (1.5*C_0), u_mean)
end

function PSP_motion_params(omega_bar::T, omega_sigma_2::T, C_0::T, B_format::String, c_phi::T, c_t::T,u_mean::T) where T<:AbstractFloat
    return PSP_params(omega_bar, omega_sigma_2, C_0, B_format, c_phi, c_t), motion_params(omega_bar,C_0, B_format, u_mean)
end

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

function set_phi_as_ic_up1!(phi_array::Array{Float64})
    #Initial_condition == "Uniform phi1"
    nparticles = size(phi_array)[2]
    phi_array[2,:] = abs.(phi_eps*randn(nparticles)) #pdf can't find zeros
    phi_array[1,:] .= 1;
end
function set_phi_as_ic_td!(phi_array::Array{Float64})
    #Initial_condition == "triple_delta"
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
function set_phi_as_ic_2l!(phi_array::Array{TF},yp::Array{TF},space_cells::CellGrid{TF}) where TF<:AbstractFloat
    #Initial_condition == "2 layers"
    nparticles = size(phi_array)[2]
    local noise_term = randn(nparticles)
    # local uniform_noise = rand(nparticles).-0.5

    phi_array[2,[yp.>0.5*height_domain]] = abs.(phi_eps*noise_term[yp.>0.5*space_cells.height_domain] )
    phi_array[1,[yp.>0.5*height_domain]] .= 1

    phi_array[1,[yp.<=0.5*height_domain]] = abs.(phi_eps*noise_term[yp.<=0.5*space_cells.height_domain] )
    phi_array[2,[yp.<=0.5*height_domain]] .= 1 #.+ uniform_noise[yp[particles,1].<=0.5*height_domain].*0.05
end
function set_phi_as_ic_dd!(phi_array::Array{Float64})
    #Initial_condition == "double delta"
    nparticles = size(phi_array)[2]
    local noise_term = randn(nparticles)
    local delta_selector = rand([true,false], nparticles)
    # local uniform_noise = rand(nparticles).-0.5

    phi_array[2,delta_selector] = abs.(phi_eps*noise_term[delta_selector] )
    phi_array[1,delta_selector] .= 1 #.+ uniform_noise[delta_selector.==1].*0.05

    phi_array[1,delta_selector] = abs.(phi_eps*noise_term[.!delta_selector] )
    phi_array[2,delta_selector] .= 1 #.+ uniform_noise[delta_selector.==2].*0.05
end

function assign_f_phi_cell!(f_phi_cell::AbstractArray{TF},phi_array::AbstractArray{TF}, psi_mesh::PsiGrid{TF}) where TF <: AbstractFloat
    "use if cell_points has alreday been determined"
    for psi1_i=1:psi_mesh.psi_partions_num
        in_1 = psi_mesh.psi_1[psi1_i].<=phi_array[1,:].<psi_mesh.psi_1[psi1_i+1]
        for psi2_i=1:psi_mesh.psi_partions_num
            in_2 = psi_mesh.psi_2[psi2_i].<=phi_array[2,:].<psi_mesh.psi_2[psi2_i+1]
            f_phi_cell[:,:] .= sum(in_1.&in_2)
        end
    end
    f_phi_cell[:, :] = f_phi_cell./(size(phi_array)[2])
end 

function assign_f_phi!(f_phi_t::Array{TF},phi_array::Array{TF}, xp::Array{TF}, yp::Array{TF}, psi_mesh::PsiGrid{TF}, space_cells::CellGrid{TF}) where TF <: AbstractFloat
    for i in 1:space_cells.y_res
        in_y = space_cells.y_edges[i].<=yp[:].<space_cells.y_edges[i+1]
        for j in 1:space_cells.x_res
            in_x = space_cells.x_edges[j].<=xp[:].<space_cells.x_edges[j+1]
            cell_particles = findall(in_x.&in_y)
            assign_f_phi_cell!(f_phi_t[:,:,i, j],phi_array[:,cell_particles],psi_mesh)
        end 
    end
end 

function bc_absorbtion_clt!(abs_points::Array{TF}, nvpart_per_part::Int, bc_k::TF, turb_k_e::Array{TF}, C_0::TF, B::TF) where TF<:AbstractFloat
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

#using binomal noise for small numbers of vparticles
function bc_absorbtion_bino!(abs_points::Array{TF}, nvpart_per_part::Int, bc_k::TF, turb_k_e::Array{TF}, C_0::TF, B::TF) where TF<:AbstractFloat
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

function bc_absorbtion_inf!(abs_points::Array{TF}, nvpart_per_part::TF, bc_k::TF, turb_k_e::Array{TF}, C_0::TF, B::TF) where TF<:AbstractFloat
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

function bc_absorbtion_bino!(abs_points::Array{TF}, nvpart_per_part::TF, bc_k::TF, turb_k_e::Array{TF}, C_0::TF, B::TF) where TF<:AbstractFloat
    bc_absorbtion_inf!(abs_points, nvpart_per_part, bc_k, turb_k_e, C_0, B) #basically a function pointer
end
function bc_absorbtion_clt!(abs_points::Array{TF}, nvpart_per_part::TF, bc_k::TF, turb_k_e::Array{TF}, C_0::TF, B::TF) where TF<:AbstractFloat
    bc_absorbtion_inf!(abs_points, nvpart_per_part, bc_k, turb_k_e, C_0, B) #basically a function pointer
end

function particle_motion_model(x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::Array{T,2}, m_params::MotionParams{T}, dt::T, space_cells::CellGrid{T}) where T<:AbstractFloat
    omega_bar=m_params.omega_bar
    C_0=m_params.C_0
    B=m_params.B
    u_mean=m_params.u_mean
    np = size(x_pos)[1]
    nt = size(x_pos)[2]-1
    bc_interact = zeros(Bool, np, nt, 4)
    #intitial vaules of velocity, maintaining consitancy with energy
    uxp = randn(np).*sqrt.(2/3 .*turb_k_e[:,1])
    uyp = randn(np).*sqrt.(2/3 .*turb_k_e[:,1])
    for t=1:nt
        for i in 1:space_cells.y_res
            in_y = space_cells.y_edges[i].<y_pos[:,t].<space_cells.y_edges[i+1]
            for j in 1:x_res
                in_x = space_cells.x_edges[j].<x_pos[:,t].<space_cells.x_edges[j+1]
                cell_particles = findall(in_x.&in_y)
                turb_k_e[cell_particles,t].=0.5*(st.mean(uxp[cell_particles].^2)+st.mean(uyp[cell_particles].^2))*1.5 #turb_e_init;
            end
        end
        print(size(uyp),' ',size(uxp),' ',size(u_mean))
        uxp = uxp+(-0.5*B*omega_bar*uxp)*dt+randn(np).*sqrt.(C_0.*turb_k_e[:,t].*omega_bar.*dt); 
        uyp = uyp+(-0.5*B*omega_bar*uyp)*dt+randn(np).*sqrt.(C_0.*turb_k_e[:,t].*omega_bar.*dt); 
        uxp_full = u_mean .+ uxp
        x_pos[:,t+1] = x_pos[:,t] + (uxp_full)*dt # random walk in x-direction
        y_pos[:,t+1] = y_pos[:,t] + uyp*dt # random walk in y-direction

            # Reflection particles at boundaries

        # Reflection at upper boundary y>height_domain
        # doing closed on top open on bottom, as cell detection is open on top,
        # closed on bottom
        mag = findall(y_pos[:,t+1].>=space_cells.height_domain) # index of particle with yp>height_domain
        dim_mag = size(mag) # dimension of array "mag"

        y_mag_succ = y_pos[mag,t+1] # yp at time t+1 corresponding to the index "mag"

        V1 = space_cells.height_domain.*ones(dim_mag) 

        ypr_mag = V1*2 .- y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
        uyp[mag] = -uyp[mag] #reflecting velocity
        bc_interact[mag,t,1] .= true

        # Reflection at lower boundary y<0
        mag = findall(y_pos[:,t+1].<=0) # index of particle with yp>height_domain
        dim_mag = size(mag) # dimension of array "mag"

        y_mag_succ = y_pos[mag,t+1] # yp at time t+1 corresponding to the index "mag"

        ypr_mag = - y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
        uyp[mag] = -uyp[mag] #reflecting velocity
        bc_interact[mag,t,2] .= true

        #bc at end (y=length_domain) of domain
        end_indicies = x_pos[:,t+1].>=space_cells.length_domain #index of particle with xp>length

        end_x = x_pos[end_indicies,t+1]
        xpr_end = end_x .- space_cells.length_domain #shifting particles back to begining
        x_pos[end_indicies,t+1] = xpr_end #replacing x coords

        bc_interact[end_indicies,t,3] .= true

        #bc at start (x=0) of domain
        start_indicies = x_pos[:,t+1].<=0 #index of particle with xp>length

        xpr_start = space_cells.length_domain .+ x_pos[start_indicies,t+1] 
        x_pos[start_indicies,t+1] = xpr_start #replacing x coords
        bc_interact[start_indicies,t,4] .= true
    end
    return bc_interact
end

function particle_motion_model(x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::T, m_params::MotionParams{T}, dt::T, space_cells::CellGrid{T}) where T<:AbstractFloat
    #for constant kinetic energy
    omega_bar=m_params.omega_bar
    C_0=m_params.C_0
    B=m_params.B
    u_mean=m_params.u_mean
    np = size(x_pos)[1]
    nt = size(x_pos)[2]-1
    bc_interact = zeros(Bool, np, nt, 4)
    #intitial vaules of velocity, maintaining consitancy with energy
    uxp = randn(np).*sqrt.(2/3 .*turb_k_e)
    uyp = randn(np).*sqrt.(2/3 .*turb_k_e)
    for t=1:nt
        print(size(uyp),' ',size(uxp),' ',size(u_mean))
        uxp = uxp+(-0.5*B*omega_bar*uxp)*dt+randn(np).*sqrt.(C_0.*turb_k_e.*omega_bar.*dt); 
        uyp = uyp+(-0.5*B*omega_bar*uyp)*dt+randn(np).*sqrt.(C_0.*turb_k_e.*omega_bar.*dt); 
        uxp_full = u_mean .+ uxp
        x_pos[:,t+1] = x_pos[:,t] + (uxp_full)*dt # random walk in x-direction
        y_pos[:,t+1] = y_pos[:,t] + uyp*dt # random walk in y-direction

            # Reflection particles at boundaries

        # Reflection at upper boundary y>height_domain
        # doing closed on top open on bottom, as cell detection is open on top,
        # closed on bottom
        mag = findall(y_pos[:,t+1].>=space_cells.height_domain) # index of particle with yp>height_domain
        dim_mag = size(mag) # dimension of array "mag"

        y_mag_succ = y_pos[mag,t+1] # yp at time t+1 corresponding to the index "mag"

        V1 = space_cells.height_domain.*ones(dim_mag) 

        ypr_mag = V1*2 .- y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
        uyp[mag] = -uyp[mag] #reflecting velocity
        bc_interact[mag,t,1] .= true

        # Reflection at lower boundary y<0
        mag = findall(y_pos[:,t+1].<=0) # index of particle with yp>height_domain
        dim_mag = size(mag) # dimension of array "mag"

        y_mag_succ = y_pos[mag,t+1] # yp at time t+1 corresponding to the index "mag"

        ypr_mag = - y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
        uyp[mag] = -uyp[mag] #reflecting velocity
        bc_interact[mag,t,2] .= true

        #bc at end (y=length_domain) of domain
        end_indicies = x_pos[:,t+1].>=space_cells.length_domain #index of particle with xp>length

        end_x = x_pos[end_indicies,t+1]
        xpr_end = end_x .- space_cells.length_domain #shifting particles back to begining
        x_pos[end_indicies,t+1] = xpr_end #replacing x coords

        bc_interact[end_indicies,t,3] .= true

        #bc at start (x=0) of domain
        start_indicies = x_pos[:,t+1].<=0 #index of particle with xp>length

        xpr_start = space_cells.length_domain .+ x_pos[start_indicies,t+1] 
        x_pos[start_indicies,t+1] = xpr_start #replacing x coords
        bc_interact[start_indicies,t,4] .= true
    end
    return bc_interact
end

function PSP_model!(f_phi::Array{TF,5},x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::Array{T,2}, bc_interact::Array{Bool,3}, bc_k::TF, dt::T, Initial_condition::String,  p_params::PSPParams{T}, psi_mesh::PsiGrid{T}, space_cells::CellGrid{T}) where T<:AbstractFloat
    omega_mean=p_params.omega_bar
    omega_sigma_2 = p_params.omega_sigma_2
    T_omega = p_params.T_omega
    C_0=p_params.C_0
    B=p_params.B
    c_phi = p_params.c_phi
    c_t = p_params.c_t
    np, nt = size(x_pos)

    phip = zeros((2, np, nt+1)) #scalar concentration at these points
    phi_pm = zeros(Int, 2, np) #pm pairs for each particle
    for t in 1:nt    
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
        for i in 1:space_cells.y_res
            in_y = space_cells.y_edges[i].<=y_pos[:,t+1].<space_cells.y_edges[i+1]
            for j in 1:space_cells.x_res
                in_x = space_cells.x_edges[j].<=x_pos[:,t+1].<space_cells.x_edges[j+1]
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
        T_mat=zeros(2,2)
        T_mat[:,1] = e_1  #coord transform matrix
        T_mat[:,2] = e_2
        dphi = T_mat\dphi  #transform to new coords
        #performing adjustment to mean 0
        corr_factor = zeros(2,np)
        for i in 1:space_cells.y_res
            in_y = space_cells.y_edges[i].<=y_pos[:,t+1].<space_cells.y_edges[i+1]
            for j in 1:space_cells.x_res
                in_x = space_cells.x_edges[j].<=x_pos[:,t+1].<space_cells.x_edges[j+1]
                cell_particles = findall(in_x.&in_y)
                for phi_i=1:2
                    phi_mean = mean(dphi[phi_i,cell_particles])
                    if phi_mean != 0 #isn't true for empty cells
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
        dphi = T_mat*dphi #return to old coords
        phip[:,:,t+1] = phip[:,:,t]+dphi
        if !(Initial_condition == "triple delta")
            phip[:,:,t+1] = phip[:,:,t+1].*(phip[:,:,t+1].>0) #forcing positive concentration
        end

        assign_f_phi!(f_phi[:,:,:,:,t],phip[:,:,t], x_pos[:,t], y_pos[:,t], psi_mesh, space_cells)
    end
end