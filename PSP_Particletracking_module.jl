module PSPParticleTrackingReactions

using Random, Distributions
import StatsBase as sb
import LinearAlgebra as la
import Statistics as st

phi_eps=1e-8
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
    c_phi::T
    c_t::T
end

struct MotionParams{T<:AbstractFloat}
    omega_bar::T
    C_0::T
    B::T
    u_mean::T
end

struct BCParams{T<:AbstractFloat, Tvp<:Real, bc_CLT}
    bc_k::T
    C_0::T
    B::T
    num_vp::Tvp
end

function cell_grid(x_res ::Int,y_res ::Int,length_domain ::T,height_domain ::T) where T<:AbstractFloat #constructor for CellGrid
    return CellGrid(x_res,y_res,length_domain,height_domain,LinRange(0,length_domain,x_res+1),LinRange(0,height_domain,y_res+1))
end

function psi_grid(psi_partions_num ::Int, phi_domain ::AbstractVector{T}) where T<:AbstractFloat #constructor for PsiGrid
    return PsiGrid(psi_partions_num, phi_domain,LinRange(phi_domain[1], phi_domain[2], psi_partions_num+1),LinRange(phi_domain[1], phi_domain[2], psi_partions_num+1))
end

function PSP_params(omega_bar::T, omega_sigma_2::T, c_phi::T, c_t::T) where T<:AbstractFloat
    return PSPParams(omega_bar, omega_sigma_2, 1/omega_bar, c_phi, c_t)
end

function motion_params(omega_bar::T,C_0::T, B_format::String, u_mean::T) where T<:AbstractFloat
    if B_format == "Decay"
        return MotionParams(omega_bar, C_0, (1+1.5*C_0), u_mean)
    elseif B_format == "Constant"
        return MotionParams(omega_bar, C_0, (1.5*C_0), u_mean)
    end
end

function BC_params(bc_k::T, C_0::T, B_format::String, num_vp::Tvp, bc_CLT::Bool) where T<:AbstractFloat where Tvp<:Int
    if B_format == "Decay"
        return BCParams{T,Tvp,bc_CLT}(bc_k, C_0, (1+1.5*C_0),num_vp)
    elseif B_format == "Constant"
        return BCParams{T,Tvp,bc_CLT}(bc_k, C_0, (1.5*C_0),num_vp)
    end
end

function BC_params(bc_k::T, C_0::T, B_format::String, num_vp::Tvp, bc_CLT::Bool) where T<:AbstractFloat where Tvp<:AbstractFloat
    num_vp == Inf || throw(DomainError("nvpart_per_part must be Int or Inf"))
    if B_format == "Decay"
        return BCParams{T,Tvp,bc_CLT}(bc_k, C_0, (1+1.5*C_0),num_vp)
    elseif B_format == "Constant"
        return BCParams{T,Tvp,bc_CLT}(bc_k, C_0, (1.5*C_0),num_vp)
    end
end

function PSP_motion_bc_params(omega_bar::T, omega_sigma_2::T, C_0::T, B_format::String, c_phi::T, c_t::T,u_mean::T,bc_k::T,num_vp::Real,bc_CLT::Bool) where T<:AbstractFloat
    return PSP_params(omega_bar, omega_sigma_2, c_phi, c_t), motion_params(omega_bar,C_0, B_format, u_mean), BC_params(bc_k, C_0, B_format, num_vp, bc_CLT)
end

function assign_pm!(phi_pm::Matrix{Int}, phi_array_t::Array{T}, particles::Vector{Int}, cell_points::Vector{Int}) where T<:AbstractFloat
    # might be worth a strategy that removes tested pairs, can't see how to make it not require a biiiig temp variable though 
    for particle in particles
        while true
            try_pm = sb.sample(cell_points,2,replace=false, ordered=true) #ordering added to reduce sample space as order for pm pairs is irrlevent
            test_dot = la.dot((phi_array_t[:,try_pm[1]]-phi_array_t[:,particle]),(phi_array_t[:,try_pm[2]]-phi_array_t[:,particle]))
            if test_dot<= 0
                phi_pm[:,particle] = try_pm
                break
            end
        end
    end
    return nothing
end

function assign_pm_single!(phi_pm::Matrix{Int}, phi_array_t::Array{T}, particles::Vector{Int}, cell_points::Vector{Int}, index_tbc::Int) where T<:AbstractFloat
    known_index = 1-(index_tbc-1)+1 #flips 1 to 2 and 2 to 1
    cell_points_set = Set(cell_points)
    for particle in particles
        !in(particle, cell_points) && throw(ArgumentError("Particles must be in cell"))
        cell_points_set_ = copy(cell_points_set)
        while true
            try_pm = rand(cell_points_set_) 
            test_dot = la.dot((phi_array_t[:,try_pm[1]]-phi_array_t[:,particle]),(phi_array_t[:,phi_pm[known_index,particle]]-phi_array_t[:,particle]))
            if  test_dot<= 0
                phi_pm[index_tbc,particle] = try_pm[1]
                break
            end
            setdiff!(cell_points_set_,try_pm)
        end
    end
    return nothing
end

function set_phi_as_ic_up1!(phi_array::Array{TF,3}, t_index::Int) where TF<:AbstractFloat
    #Initial_condition == "Uniform phi1"
    nparticles = size(phi_array)[2]
    phi_array[2,:,t_index] = abs.(phi_eps*randn(TF, nparticles)) #pdf can't find zeros
    phi_array[1,:,t_index] .= 1 
    return nothing
end
#using subset of particles, mainly for inflow
function set_phi_as_ic_up1!(phi_array::Array{TF,3}, t_index::Int, subset_indicies) where TF<:AbstractFloat
    #Initial_condition == "Uniform phi1"
    nparticles = size(phi_array[:,subset_indicies,:])[2]
    phi_array[2,subset_indicies,t_index] = abs.(phi_eps*randn(TF, nparticles)) #pdf can't find zeros
    phi_array[1,subset_indicies,t_index] .= 1 
    return nothing
end
function set_phi_as_ic_td!(phi_array::Array{TF,3}, t_index::Int) where TF<:AbstractFloat
    #Initial_condition == "triple_delta"
    nparticles = size(phi_array)[2]
    local delta_selector = rand(1:3, nparticles)
    local noise_term = randn(TF, nparticles)

    phi_array[1,delta_selector.=1,t_index] = -sqrt(3)/2 .+phi_eps .*noise_term[delta_selector.=1]
    phi_array[2,delta_selector.=1,t_index] .= -0.5

    phi_array[1,delta_selector.=2,t_index] = sqrt(3)/2 .+phi_eps .*noise_term[delta_selector.=2]
    phi_array[2,delta_selector.=2,t_index] .= -0.5

    phi_array[1,delta_selector.=3,t_index] .= phi_eps
    phi_array[2,delta_selector.=3,t_index] = 1.0 .+phi_eps .*noise_term[delta_selector.=3]
    return nothing
end
function set_phi_as_ic_2l!(phi_array::Array{TF,3},yp::Vector{TF},space_cells::CellGrid{TF}, t_index::Int) where TF<:AbstractFloat
    #Initial_condition == "2 layers"
    nparticles = size(phi_array)[2]
    local noise_term = randn(TF, nparticles)
    # local uniform_noise = rand(nparticles).-0.5

    phi_array[2,[yp.>0.5*height_domain],t_index] = abs.(phi_eps*noise_term[yp.>0.5*space_cells.height_domain] )
    phi_array[1,[yp.>0.5*height_domain],t_index] .= 1

    phi_array[1,[yp.<=0.5*height_domain],t_index] = abs.(phi_eps*noise_term[yp.<=0.5*space_cells.height_domain] )
    phi_array[2,[yp.<=0.5*height_domain],t_index] .= 1 #.+ uniform_noise[yp[particles,1].<=0.5*height_domain].*0.05
    return nothing
end
function set_phi_as_ic_dd!(phi_array::Array{TF,3},t_index::Int) where TF<:AbstractFloat
    #Initial_condition == "double delta"
    nparticles = size(phi_array)[2]
    local noise_term = randn(TF, nparticles)
    local delta_selector = rand([true,false], nparticles)
    # local uniform_noise = rand(nparticles).-0.5

    phi_array[2,delta_selector,t_index] = abs.(phi_eps*noise_term[delta_selector] )
    phi_array[1,delta_selector,t_index] .= 1 #.+ uniform_noise[delta_selector.==1].*0.05

    phi_array[1,delta_selector,t_index] = abs.(phi_eps*noise_term[.!delta_selector] )
    phi_array[2,delta_selector,t_index] .= 1 #.+ uniform_noise[delta_selector.==2].*0.05
    return nothing
end

function assign_f_phi_cell!(f_phi_cell::AbstractArray{TF,5},phi_array::AbstractArray{TF,2}, psi_mesh::PsiGrid{TF}, cell_row::Int, cell_column::Int, t_index::Int) where TF <: AbstractFloat
    "use if cell_points has alreday been determined"
    phi_array_sort_ind_1 = sortperm(phi_array[1,:])
    phi1_i = 2 #upper bound index #assuming no values outside range
    prev_stop_1 = 1
    current_i = 0
    for i in phi_array_sort_ind_1
        current_i += 1
        while phi_array[1,i] > psi_mesh.psi_1[phi1_i]
            phi_array_sort_ind_2 = phi_array_sort_ind_1[sortperm(phi_array[2,phi_array_sort_ind_1[prev_stop_1:current_i-1]]).+prev_stop_1.-1]
            accum = 0
            phi2_i = 2
            for j in phi_array_sort_ind_2
                while phi_array[2,j] > psi_mesh.psi_2[phi2_i]
                    f_phi_cell[phi1_i-1, phi2_i-1, cell_row, cell_column, t_index] = float(accum)
                    accum = 0
                    phi2_i += 1
                end
                accum +=1
            end
            try 
                f_phi_cell[phi1_i-1, phi2_i-1, cell_row, cell_column, t_index] = float(accum)
            catch 
                println(phi1_i,' ', phi2_i,' ', cell_row,' ', cell_column,' ', t_index)
                println(maximum(phi_array[1,:]),' ',maximum(phi_array[2,:]))
                println(psi_mesh.phi_domain)
                rethrow()
            end
            prev_stop_1 = current_i
            phi1_i +=1
        end
    end
    #catch those for which there are no higher points
    phi_array_sort_ind_2 = phi_array_sort_ind_1[sortperm(phi_array[2,phi_array_sort_ind_1[prev_stop_1:current_i]]).+prev_stop_1.-1]
    accum = 0
    phi2_i = 2
    for j in phi_array_sort_ind_2
        while phi_array[2,j] > psi_mesh.psi_2[phi2_i]
            f_phi_cell[phi1_i-1, phi2_i-1, cell_row, cell_column, t_index] = float(accum)
            accum = 0
            phi2_i += 1
        end
        accum +=1
    end
    f_phi_cell[phi1_i-1, phi2_i-1, cell_row, cell_column, t_index] = float(accum)

    f_phi_cell[:, :,cell_row, cell_column, t_index] = f_phi_cell[:,:,cell_row, cell_column, t_index]./(size(phi_array)[2])
    any(f_phi_cell[:, :, cell_row, cell_column, t_index] .!= f_phi_cell[:, :, cell_row, cell_column, t_index]) && (print("error",' ',size(phi_array),' '), qwertuio)
    return nothing
end 

function assign_f_phi!(f_phi_t::Array{TF},phi_array::Array{TF}, xp::Vector{TF}, yp::Vector{TF}, psi_mesh::PsiGrid{TF}, space_cells::CellGrid{TF}, t_index::Int) where TF <: AbstractFloat
    x_sort = sortperm(xp)
    x_i = 2 #assuming no values outside range
    prev_stop_x = 1
    current_i = 0
    for i in x_sort
        current_i += 1
        while xp[i] > space_cells.x_edges[x_i]
            y_sort = x_sort[sortperm(yp[x_sort[prev_stop_x:current_i-1]]).+(prev_stop_x-1)]
            prev_stop_y = 1
            y_i = 2
            current_j = 0
            for j in y_sort
                current_j += 1
                while yp[j] > space_cells.y_edges[y_i]
                    assign_f_phi_cell!(f_phi_t, phi_array[:,y_sort[prev_stop_y:current_j-1]], psi_mesh, y_i-1, x_i-1, t_index)
                    prev_stop_y = current_j
                    y_i+=1
                end
            end
            assign_f_phi_cell!(f_phi_t, phi_array[:,y_sort[prev_stop_y:current_j]], psi_mesh, y_i-1, x_i-1, t_index)
            prev_stop_x = current_i
            x_i +=1
        end
    end
    #catch final cell
    y_sort = x_sort[sortperm(yp[x_sort[prev_stop_x:current_i]]).+(prev_stop_x-1)]
    prev_stop_y = 1
    y_i = 2
    current_j = 0
    for j in y_sort
        current_j += 1
        while yp[j] > space_cells.y_edges[y_i]
            assign_f_phi_cell!(f_phi_t, phi_array[:,y_sort[prev_stop_y:current_j-1]], psi_mesh, y_i-1, x_i-1, t_index)
            prev_stop_y = current_j
            y_i+=1
        end
    end
    assign_f_phi_cell!(f_phi_t, phi_array[:,y_sort[prev_stop_y:current_j]], psi_mesh, y_i-1, x_i-1, t_index)
    return nothing
end 

function eval_by_cell!(func!::Function, xp::Vector{TF}, yp::Vector{TF}, space_cells::CellGrid{TF}) where TF <: AbstractFloat
    "func! is expected to be a function of signature func(cell_row,cell_column,cell_particles)->nothing"
    x_sort = sortperm(xp)
    x_i = 2 #assuming no values outside range
    prev_stop_x = 1
    current_i = 0
    for i in x_sort
        current_i += 1
        while xp[i] > space_cells.x_edges[x_i]
            y_sort = x_sort[sortperm(yp[x_sort[prev_stop_x:current_i-1]]).+(prev_stop_x-1)]
            prev_stop_y = 1
            y_i = 2
            current_j = 0
            for j in y_sort
                current_j += 1
                while yp[j] > space_cells.y_edges[y_i]
                    func!(y_i-1,x_i-1,y_sort[prev_stop_y:current_j-1])#calling function
                    prev_stop_y = current_j
                    y_i+=1
                end
            end
            func!(y_i-1,x_i-1,y_sort[prev_stop_y:current_j])#calling function
            prev_stop_x = current_i
            x_i +=1
        end
    end
    #catch final cell
    y_sort = x_sort[sortperm(yp[x_sort[prev_stop_x:current_i]]).+(prev_stop_x-1)]
    prev_stop_y = 1
    y_i = 2
    current_j = 0
    for j in y_sort
        current_j += 1
        while yp[j] > space_cells.y_edges[y_i]
            func!(y_i-1,x_i-1,y_sort[prev_stop_y:current_j-1])#calling function
            prev_stop_y = current_j
            y_i+=1
        end
    end
    func!(y_i-1,x_i-1,y_sort[prev_stop_y:current_j])#calling function
    return nothing
end 

#CLt/normal
function bc_absorbtion!(phip::Array{TF,3}, abs_points::Vector{Bool}, turb_k_e::Vector{TF}, bc_params::BCParams{TF,Int,true}, t_index::Int) where TF<:AbstractFloat
    n_abs = sum(abs_points)
    abs_k = bc_params.bc_k.*ones(TF, 2,n_abs)
    effective_v_particles =( phip[:,abs_points,t_index].*bc_params.num_vp)
    #K for Erban and Chapman approximation 
    P = zeros(TF, 2,n_abs)
    P[1,:] = min.(abs_k[1,:].*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)
    P[2,:] = min.(abs_k[2,:].*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)
    #by CLT approx dist for number of virtual particles to have reacted
    xi = randn(TF, 2,n_abs).*sqrt.((P.*(1 .-P)))
    #catching places where all mass has been removed
    xi = [effective_v_particles[i,j]>0 ? xi[i,j]./sqrt(effective_v_particles[i,j]) : 0 for i in 1:2, j in 1:n_abs]
    ratios = max.(min.((1 .-P) + xi,1),0)
    phip[:, abs_points, t_index] = phip[:, abs_points, t_index].*ratios
    return nothing
end

#CLt/normal Precomp
function bc_absorbtion!(phip::Array{TF,3}, abs_points::Vector{Bool}, bc_params::BCParams{TF,Int,true}, t_index::Int, Precomp_P::TF) where TF<:AbstractFloat
    n_abs = sum(abs_points)
    effective_v_particles =( phip[:,abs_points,t_index].*bc_params.num_vp)
    #K for Erban and Chapman approximation 
    P = zeros(TF, 2,n_abs)
    P .= Precomp_P
    #by CLT approx dist for number of virtual particles to have reacted
    xi = randn(TF, 2,n_abs).*sqrt.((P.*(1 .-P)))
    #catching places where all mass has been removed
    xi = [effective_v_particles[i,j]>0 ? xi[i,j]./sqrt(effective_v_particles[i,j]) : 0 for i in 1:2, j in 1:n_abs]
    ratios = max.(min.((1 .-P) + xi,1),0)
    phip[:, abs_points, t_index] = phip[:, abs_points, t_index].*ratios
    return nothing
end

#using binomal noise for small numbers of vparticles
function bc_absorbtion!(phip::Array{TF,3}, abs_points::Vector{Bool}, turb_k_e::Vector{TF}, bc_params::BCParams{TF,Int,false}, t_index::Int) where TF<:AbstractFloat
    n_abs = sum(abs_points)
    abs_k = bc_params.bc_k.*ones(TF,2,n_abs)
    effective_v_particles =( phip[:,abs_points,t_index].*bc_params.num_vp)
    #K for Erban and Chapman approximation 
    P = zeros(TF,2,n_abs)
    P[1,:] = min.(abs_k[1,:].*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)
    P[2,:] = min.(abs_k[2,:].*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)
    #Binomal dist for number of virtual particles to have reacted
    xi_dist = Binomial.(TF, ceil.(effective_v_particles),1 .-P)
    xi = [rand(xi_dist[i,j]) for i in 1:2, j in 1:n_abs]
    ratios = [effective_v_particles[i,j]>0 ? xi[i,j]./ceil.(effective_v_particles[i,j]) : 0 for i in 1:2, j in 1:n_abs]
    phip[:, abs_points, t_index] = phip[:, abs_points, t_index].*ratios
    return nothing
end

#binomal prcomp
function bc_absorbtion!(phip::Array{TF,3}, abs_points::Vector{Bool}, bc_params::BCParams{TF,Int,false}, t_index::Int, Precomp_P::TF) where TF<:AbstractFloat
    n_abs = sum(abs_points)
    effective_v_particles =( phip[:,abs_points,t_index].*bc_params.num_vp)
    #K for Erban and Chapman approximation 
    P = zeros(TF,2,n_abs)
    P .= Precomp_P
    #Binomal dist for number of virtual particles to have reacted
    xi_dist = Binomial.(TF,ceil.(effective_v_particles),1 .-P)
    xi = [rand(xi_dist[i,j]) for i in 1:2, j in 1:n_abs]
    ratios = [effective_v_particles[i,j]>0 ? xi[i,j]./ceil.(effective_v_particles[i,j]) : 0 for i in 1:2, j in 1:n_abs]
    phip[:, abs_points, t_index] = phip[:, abs_points, t_index].*ratios
    return nothing
end

#mean
function bc_absorbtion!(phip::Array{TF,3}, abs_points::Vector{Bool}, turb_k_e::Vector{TF}, bc_params::BCParams{TF,TF}, t_index::Int) where TF<:AbstractFloat 
    n_abs = sum(abs_points)
    abs_k = bc_params.bc_k.*ones(2,n_abs)
    #K for Erban and Chapman approximation 
    P = zeros(TF,2,n_abs)
    P[1,:] = min.(abs_k[1,:].*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)
    P[2,:] = min.(abs_k[2,:].*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)
    ratios = 1 .-P #taking mean for limiting case
    phip[:, abs_points, t_index] = phip[:, abs_points, t_index].*ratios
    return nothing
end

#mean Precomp
function bc_absorbtion!(phip::Array{TF,3}, abs_points::Vector{Bool}, bc_params::BCParams{TF,TF}, t_index::Int, Precomp_P::TF) where TF<:AbstractFloat 
    #K for Erban and Chapman approximation 
    ratios = 1 - Precomp_P #taking mean for limiting case
    phip[:, abs_points, t_index] = phip[:, abs_points, t_index].*ratios
    return nothing
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
    uxp = randn(T, np).*sqrt.(2/3 .*turb_k_e[:,1])
    uyp = randn(T, np).*sqrt.(2/3 .*turb_k_e[:,1])
    for t=1:nt
        uxp = uxp+(-0.5*B*omega_bar*uxp)*dt+randn(T, np).*sqrt.(C_0.*turb_k_e[:,t].*omega_bar.*dt); 
        uyp = uyp+(-0.5*B*omega_bar*uyp)*dt+randn(T, np).*sqrt.(C_0.*turb_k_e[:,t].*omega_bar.*dt); 
        for i in 1:space_cells.y_res
            in_y = space_cells.y_edges[i].<y_pos[:,t].<space_cells.y_edges[i+1]
            for j in 1:x_res
                in_x = space_cells.x_edges[j].<x_pos[:,t].<space_cells.x_edges[j+1]
                cell_particles = findall(in_x.&in_y)
                turb_k_e[cell_particles,t+1].=0.5*(st.mean(uxp[cell_particles].^2)+st.mean(uyp[cell_particles].^2))*1.5 #turb_e_init;
            end
        end
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

        V1 = space_cells.height_domain.*ones(T, dim_mag) 

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

#there's a burn in here now
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
    uxp = randn(T, np).*sqrt.(2/3 .*turb_k_e)
    uyp = randn(T, np).*sqrt.(2/3 .*turb_k_e)
    burn_in=10

    #ensures the velocities have the correct distribution 
    for t=1:burn_in
        uxp = uxp+(-0.5*B*omega_bar*uxp)*dt+randn(T, np).*sqrt.(C_0.*turb_k_e.*omega_bar.*dt); 
        uyp = uyp+(-0.5*B*omega_bar*uyp)*dt+randn(T, np).*sqrt.(C_0.*turb_k_e.*omega_bar.*dt); 
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

        V1 = space_cells.height_domain.*ones(T, dim_mag) 

        ypr_mag = V1*2 .- y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
        uyp[mag] = -uyp[mag] #reflecting velocity

        # Reflection at lower boundary y<0
        mag = findall(y_pos[:,t+1].<=0) # index of particle with yp>height_domain
        dim_mag = size(mag) # dimension of array "mag"

        y_mag_succ = y_pos[mag,t+1] # yp at time t+1 corresponding to the index "mag"

        ypr_mag = - y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp<0 with yp of reflected particle
        uyp[mag] = -uyp[mag] #reflecting velocity

        #bc at end (y=length_domain) of domain
        end_indicies = x_pos[:,t+1].>=space_cells.length_domain #index of particle with xp>length

        end_x = x_pos[end_indicies,t+1]
        xpr_end = end_x .- space_cells.length_domain #shifting particles back to begining
        x_pos[end_indicies,t+1] = xpr_end #replacing x coords


        #bc at start (x=0) of domain
        start_indicies = x_pos[:,t+1].<=0 #index of particle with xp>length

        xpr_start = space_cells.length_domain .+ x_pos[start_indicies,t+1] 
        x_pos[start_indicies,t+1] = xpr_start #replacing x coords
    end
    x_pos[:,1]=x_pos[:,burn_in]
    y_pos[:,1]=y_pos[:,burn_in]
    for t=1:nt
        uxp = uxp+(-0.5*B*omega_bar*uxp)*dt+randn(T, np).*sqrt.(C_0.*turb_k_e.*omega_bar.*dt); 
        uyp = uyp+(-0.5*B*omega_bar*uyp)*dt+randn(T, np).*sqrt.(C_0.*turb_k_e.*omega_bar.*dt); 
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

        V1 = space_cells.height_domain.*ones(T, dim_mag) 

        ypr_mag = V1*2 .- y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
        uyp[mag] = -uyp[mag] #reflecting velocity
        bc_interact[mag,t,1] .= true

        # Reflection at lower boundary y<0
        mag = findall(y_pos[:,t+1].<=0) # index of particle with yp>height_domain
        dim_mag = size(mag) # dimension of array "mag"

        y_mag_succ = y_pos[mag,t+1] # yp at time t+1 corresponding to the index "mag"

        ypr_mag = - y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp<0 with yp of reflected particle
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

#reflective bc at begining
function particle_motion_model_ref_start(x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::T, m_params::MotionParams{T}, dt::T, space_cells::CellGrid{T}) where T<:AbstractFloat
    #for constant kinetic energy
    omega_bar=m_params.omega_bar
    C_0=m_params.C_0
    B=m_params.B
    u_mean=m_params.u_mean
    np = size(x_pos)[1]
    nt = size(x_pos)[2]-1
    bc_interact = zeros(Bool, np, nt, 4)
    #intitial vaules of velocity, maintaining consitancy with energy
    uxp = randn(T, np).*sqrt.(2/3 .*turb_k_e)
    uyp = randn(T, np).*sqrt.(2/3 .*turb_k_e)
    for t=1:nt
        uxp = uxp+(-0.5*B*omega_bar*uxp)*dt+randn(T, np).*sqrt.(C_0.*turb_k_e.*omega_bar.*dt); 
        uyp = uyp+(-0.5*B*omega_bar*uyp)*dt+randn(T, np).*sqrt.(C_0.*turb_k_e.*omega_bar.*dt); 
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

        V1 = space_cells.height_domain.*ones(T, dim_mag) 

        ypr_mag = V1*2 .- y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
        uyp[mag] = -uyp[mag] #reflecting velocity
        bc_interact[mag,t,1] .= true

        # Reflection at lower boundary y<0
        mag = findall(y_pos[:,t+1].<=0) # index of particle with yp>height_domain
        dim_mag = size(mag) # dimension of array "mag"

        y_mag_succ = y_pos[mag,t+1] # yp at time t+1 corresponding to the index "mag"

        ypr_mag = - y_mag_succ  # yp at time t+1 of the reflected particle

        y_pos[mag,t+1]= ypr_mag #replacement of yp<0 with yp of reflected particle
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

        xpr_start = .- x_pos[start_indicies,t+1] 
        x_pos[start_indicies,t+1] = xpr_start #replacing x coords
        uxp[start_indicies] = -uxp[start_indicies] #reflecting velocity
        bc_interact[start_indicies,t,4] .= true
    end
    return bc_interact
end

#variable turb_k_e
function PSP_model!(f_phi::Array{T,5},x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::Array{T,2}, bc_interact::Array{Bool,3}, dt::T, initial_condition::String,  p_params::PSPParams{T}, psi_mesh::PsiGrid{T}, space_cells::CellGrid{T}, bc_params::BCParams{T}, verbose::Bool=false) where T<:AbstractFloat
    omega_mean=p_params.omega_bar
    omega_sigma_2 = p_params.omega_sigma_2
    T_omega = p_params.T_omega
    c_phi = p_params.c_phi
    c_t = p_params.c_t
    np, nt = size(x_pos)
    nt-=1

    phip = zeros((2, np, nt+1)) #scalar concentration at these points
    phi_pm = zeros(Int, 2, np) #pm pairs for each particle

    if initial_condition == "Uniform phi_1"
        set_phi_as_ic_up1!(phip,1)
    elseif initial_condition == "triple delta"
        set_phi_as_ic_td!(phip,1)
    elseif initial_condition == "2 layers"
        set_phi_as_ic_2l!(phip,y_pos[:,1],space_cells,1)
    elseif initial_condition == "double delta"
        set_phi_as_ic_dd!(phip,1)
    else
        throw(ArgumentError("Not a valid intitial condition"))
    end
    assign_f_phi!(f_phi,phip[:,:,1], x_pos[:,1], y_pos[:,1], psi_mesh, space_cells,1)

    omegap = zeros(T, np,nt+1) #turbulence frequency
    omega0_dist = Gamma(1/(omega_sigma_2),(omega_sigma_2)*omega_mean) #this should now match long term distribution of omega
    omegap[:,1] = rand(T, omega0_dist, np)
    
    eval_by_cell!((i,j,cell_particles)-> (assign_pm!(phi_pm, phip[:,:,1], cell_particles, cell_particles)
    ;assign_f_phi_cell!(f_phi,phip[:,cell_particles,1],psi_mesh,i,j,1);return nothing) , x_pos[:,1], y_pos[:,1], space_cells)
    # for i in 1:space_cells.y_res
    #     in_y = space_cells.y_edges[i].<y_pos[:,1].<space_cells.y_edges[i+1]
    #     for j in 1:space_cells.x_res
    #         in_x = space_cells.x_edges[j].<x_pos[:,1].<space_cells.x_edges[j+1]
    #         cell_particles = findall(in_x.&in_y)
    #         assign_pm!(phi_pm, phip[:,:,1], cell_particles, cell_particles)
    #         assign_f_phi_cell!(f_phi,phip[:,cell_particles,1],psi_mesh,i,j,1)
    #     end 
    # end

    #time stamp until new p/m bound found, needed to ensure particles are
    #decorreltaed
    t_decorr_p = 1 ./(c_t.*omegap[phi_pm[1,:],1])
    t_decorr_m = 1 ./(c_t.*omegap[phi_pm[2,:],1])

    for t in 1:nt
        verbose && print(t,' ')
        # print(maximum(phip[:,:,t]),' ')
        #E-M solver for omega 
        dw = sqrt(dt).*randn(T, np) #random draws
        omegap[:,t+1] = omegap[:,t]-(omegap[:,t].-omega_mean)./T_omega.*dt + sqrt.(omegap[:,t].*(2*omega_sigma_2*omega_mean/T_omega)).*dw
        omegap[:,t+1] = omegap[:,t+1].*(omegap[:,t+1].>0) #enforcing positivity

        #stepping the decorrelation times
        t_decorr_p = t_decorr_p.-dt;
        t_decorr_m = t_decorr_m.-dt;
        #sets index of those to be renewed with 0 - which doesn't correspond to any particle
        phi_pm[1,(t_decorr_p.<=0)] .= 0 
        phi_pm[2,(t_decorr_m.<=0)] .= 0

        #split into cells, compute centres/targets, run ODE step
        eval_by_cell!(function (i,j,cell_particles)
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
            return nothing
        end, x_pos[:,t], y_pos[:,t], space_cells)

        phi_c = 0.5.*(phip[:,phi_pm[1,:],t]+phip[:,phi_pm[2,:],t])
        diffusion = zeros(2,np)
        diffusion[1,:] = (phip[1,:,t]-phi_c[1,:]).*(exp.(-c_phi.*0.5.*omegap[:,t].*dt).-1.0)
        diffusion[2,:] = (phip[2,:,t]-phi_c[2,:]).*(exp.(-c_phi.*0.5.*omegap[:,t].*dt).-1.0)
        reaction = zeros(T, 2,np) # body reaction
        # reaction = dt.*(reaction).*exp.(c_phi.*0.5.*omegap[:,t].*dt) #integration of reation term to match diffusion scheme, uncomment if reaction !=0
        dphi = (diffusion .+ reaction)

        # ensuring mean 0 change
        # generating a random orthonormal basis
        # is 2-d so genrate a random unit vector from an angle and proceed based
        # on that
        angle = 2*pi*rand(1)[1];
        e_1 = [cos(angle),sin(angle)]
        handedness = sb.sample([-1,1],1)[1] #randomly choose betwen left or right handed system
        e_2 = handedness*[e_1[2],-e_1[1]]
        T_mat=zeros(T,2,2)
        T_mat[:,1] = e_1  #coord transform matrix
        T_mat[:,2] = e_2
        dphi = T_mat\dphi  #transform to new coords
        #performing adjustment to mean 0
        corr_factor = zeros(T, 2,np)
        
        eval_by_cell!(function (i,j,cell_particles)
            for phi_i=1:2
                phi_mean = mean(dphi[phi_i,cell_particles])
                if phi_mean != 0 #isn't true for empty cells
                    cell_points_pos = dphi[phi_i,cell_particles].>0
                    cell_points_neg = dphi[phi_i,cell_particles].<0
                    phi_pos_mean = mean(cell_points_pos.*dphi[phi_i,cell_particles])
                    phi_neg_mean = mean((cell_points_neg).*dphi[phi_i,cell_particles])
                    if phi_mean>0
                        corr_factor[phi_i,cell_particles]=.-cell_points_pos*(phi_neg_mean./phi_pos_mean) + (1 .- cell_points_pos)
                    else
                        corr_factor[phi_i,cell_particles]=.-(cell_points_neg)*(phi_pos_mean./phi_neg_mean) + (1 .- cell_points_neg)
                    end
                end
            end
            return nothing
        end, x_pos[:,t], y_pos[:,t], space_cells)
        dphi = corr_factor.*dphi
        dphi = T_mat*dphi #return to old coords
        phip[:,:,t+1] = phip[:,:,t]+dphi
        if !(initial_condition == "triple delta")
            phip[:,:,t+1] = phip[:,:,t+1].*(phip[:,:,t+1].>0) #forcing positive concentration
        end

        bc_absorbtion!(phip,bc_interact[:,t,2],turb_k_e[bc_interact[:,t,2],t+1],bc_params,t+1) #currently only reacting on bottom bc

        assign_f_phi!(f_phi,phip[:,:,t+1], x_pos[:,t+1], y_pos[:,t+1], psi_mesh, space_cells,t+1)
        # print(maximum(phip[:,:,t]),' ')
    end
    verbose && println("end")
    return nothing
end

#constant turb_k_e
function PSP_model!(f_phi::Array{T,5},x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::T, bc_interact::Array{Bool,3}, dt::T, initial_condition::String,  p_params::PSPParams{T}, psi_mesh::PsiGrid{T}, space_cells::CellGrid{T}, bc_params::BCParams{T}, verbose::Bool=false) where T<:AbstractFloat
    omega_mean=p_params.omega_bar
    omega_sigma_2 = p_params.omega_sigma_2
    T_omega = p_params.T_omega
    c_phi = p_params.c_phi
    c_t = p_params.c_t
    np, nt = size(x_pos)
    nt-=1
    precomp_P = min.(bc_params.bc_k.*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)

    phip = zeros(T, (2, np, nt+1)) #scalar concentration at these points
    phi_pm = zeros(Int, 2, np) #pm pairs for each particle

    if initial_condition == "Uniform phi_1"
        set_phi_as_ic_up1!(phip,1)
    elseif initial_condition == "triple delta"
        set_phi_as_ic_td!(phip,1)
    elseif initial_condition == "2 layers"
        set_phi_as_ic_2l!(phip,y_pos[:,1],space_cells,1)
    elseif initial_condition == "double delta"
        set_phi_as_ic_dd!(phip,1)
    else
        throw(ArgumentError("Not a valid intitial condition"))
    end
    assign_f_phi!(f_phi,phip[:,:,1], x_pos[:,1], y_pos[:,1], psi_mesh, space_cells,1)

    omegap = zeros(T, np,nt+1) #turbulence frequency
    omega0_dist = Gamma(1/(omega_sigma_2),(omega_sigma_2)*omega_mean) #this should now match long term distribution of omega
    omegap[:,1] = rand(omega0_dist, T, np)
    
    eval_by_cell!((i,j,cell_particles)-> (assign_pm!(phi_pm, phip[:,:,1], cell_particles, cell_particles)
    ;assign_f_phi_cell!(f_phi,phip[:,cell_particles,1],psi_mesh,i,j,1);return nothing) , x_pos[:,1], y_pos[:,1], space_cells)

    #time stamp until new p/m bound found, needed to ensure particles are
    #decorreltaed
    t_decorr_p = 1 ./(c_t.*omegap[phi_pm[1,:],1])
    t_decorr_m = 1 ./(c_t.*omegap[phi_pm[2,:],1])

    for t in 1:nt
        verbose && print(t,' ')
        #E-M solver for omega 
        dw = sqrt(dt).*randn(T, np) #random draws
        omegap[:,t+1] = omegap[:,t]-(omegap[:,t].-omega_mean)./T_omega.*dt + sqrt.(omegap[:,t].*(2*omega_sigma_2*omega_mean/T_omega)).*dw
        omegap[:,t+1] = omegap[:,t+1].*(omegap[:,t+1].>0) #enforcing positivity

        #stepping the decorrelation times
        t_decorr_p = t_decorr_p.-dt;
        t_decorr_m = t_decorr_m.-dt;
        #sets index of those to be renewed with 0 - which doesn't correspond to any particle
        phi_pm[1,(t_decorr_p.<=0)] .= 0 #the t_decorrs are 2d for some reason
        phi_pm[2,(t_decorr_m.<=0)] .= 0

        #split into cells, compute centres/targets, run ODE step
        eval_by_cell!(function (i,j,cell_particles)
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
            return nothing
        end, x_pos[:,t], y_pos[:,t], space_cells)

        phi_c = 0.5.*(phip[:,phi_pm[1,:],t]+phip[:,phi_pm[2,:],t])
        diffusion = zeros(T, 2,np)
        diffusion[1,:] = (phip[1,:,t]-phi_c[1,:]).*(exp.(-c_phi.*0.5.*omegap[:,t].*dt).-1.0)
        diffusion[2,:] = (phip[2,:,t]-phi_c[2,:]).*(exp.(-c_phi.*0.5.*omegap[:,t].*dt).-1.0)
        reaction = zeros(T, 2,np) # body reaction
        # reaction = dt.*(reaction).*exp.(c_phi.*0.5.*omegap[:,t].*dt) #integration of reation term to match diffusion scheme, uncomment if reaction !=0
        dphi = (diffusion .+ reaction)

        # ensuring mean 0 change
        # generating a random orthonormal basis
        # is 2-d so genrate a random unit vector from an angle and proceed based
        # on that
        angle = 2*pi*rand(T, 1)[1];
        e_1 = [cos(angle),sin(angle)]
        handedness = sb.sample([-1,1],1)[1] #randomly choose betwen left or right handed system
        e_2 = handedness*[e_1[2],-e_1[1]]
        T_mat=zeros(2,2)
        T_mat[:,1] = e_1  #coord transform matrix
        T_mat[:,2] = e_2
        dphi = T_mat\dphi  #transform to new coords
        #performing adjustment to mean 0
        corr_factor = zeros(T, 2,np)
        
        eval_by_cell!(function (i,j,cell_particles)
            for phi_i=1:2
                phi_mean = mean(dphi[phi_i,cell_particles])
                if phi_mean != 0 #isn't true for empty cells
                    cell_points_pos = dphi[phi_i,cell_particles].>0
                    cell_points_neg = dphi[phi_i,cell_particles].<0
                    phi_pos_mean = mean(cell_points_pos.*dphi[phi_i,cell_particles])
                    phi_neg_mean = mean((cell_points_neg).*dphi[phi_i,cell_particles])
                    if phi_mean>0
                        corr_factor[phi_i,cell_particles]=.-cell_points_pos*(phi_neg_mean./phi_pos_mean) + (1 .- cell_points_pos)
                    else
                        corr_factor[phi_i,cell_particles]=.-(cell_points_neg)*(phi_pos_mean./phi_neg_mean) + (1 .- cell_points_neg)
                    end
                end
            end
            return nothing
        end, x_pos[:,t], y_pos[:,t], space_cells)

        dphi = corr_factor.*dphi
        dphi = T_mat*dphi #return to old coords
        phip[:,:,t+1] = phip[:,:,t]+dphi
        if !(initial_condition == "triple delta")
            phip[:,:,t+1] = phip[:,:,t+1].*(phip[:,:,t+1].>0) #forcing positive concentration
        end

        bc_absorbtion!(phip,bc_interact[:,t,2],bc_params,t+1, precomp_P) #currently only reacting on bottom bc

        assign_f_phi!(f_phi,phip[:,:,t+1], x_pos[:,t+1], y_pos[:,t+1], psi_mesh, space_cells,t+1)
        # print(maximum(phip[:,:,t]),' ')
    end
    verbose && println("end")
    return nothing
end

function PSP_model_inflow_record_flux!(f_phi::Array{T,5},y_0_flux::Array{T,5},x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::T, bc_interact::Array{Bool,3}, dt::T, initial_condition::String,  p_params::PSPParams{T}, psi_mesh::PsiGrid{T}, space_cells::CellGrid{T}, bc_params::BCParams{T}, verbose::Bool=false) where T<:AbstractFloat
    omega_mean=p_params.omega_bar
    omega_sigma_2 = p_params.omega_sigma_2
    T_omega = p_params.T_omega
    c_phi = p_params.c_phi
    c_t = p_params.c_t
    np, nt = size(x_pos)
    nt-=1
    precomp_P = min.(bc_params.bc_k.*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)

    phip = zeros(T, (2, np, nt+1)) #scalar concentration at these points
    phi_pm = zeros(Int, 2, np) #pm pairs for each particle

    if initial_condition == "Uniform phi_1"
        set_phi_as_ic_up1!(phip,1)
    elseif initial_condition == "triple delta"
        set_phi_as_ic_td!(phip,1)
    elseif initial_condition == "2 layers"
        set_phi_as_ic_2l!(phip,y_pos[:,1],space_cells,1)
    elseif initial_condition == "double delta"
        set_phi_as_ic_dd!(phip,1)
    else
        throw(ArgumentError("Not a valid intitial condition"))
    end

    omegap = zeros(T, np,nt+1) #turbulence frequency
        omega0_dist = Gamma(1/(omega_sigma_2),(omega_sigma_2)*omega_mean) #this should now match long term distribution of omega
    omegap[:,1] = rand(omega0_dist,T, np)
    
    eval_by_cell!((i,j,cell_particles)-> (assign_pm!(phi_pm, phip[:,:,1], cell_particles, cell_particles)
    ;assign_f_phi_cell!(f_phi,phip[:,cell_particles,1],psi_mesh,i,j,1);return nothing) , x_pos[:,1], y_pos[:,1], space_cells)

    #time stamp until new p/m bound found, needed to ensure particles are
    #decorreltaed
    t_decorr_p = 1 ./(c_t.*omegap[phi_pm[1,:],1])
    t_decorr_m = 1 ./(c_t.*omegap[phi_pm[2,:],1])

    for t in 1:nt
        verbose && print(t,' ')
        #E-M solver for omega 
        dw = sqrt(dt).*randn(T, np) #random draws
        omegap[:,t+1] = omegap[:,t]-(omegap[:,t].-omega_mean)./T_omega.*dt + sqrt.(omegap[:,t].*(2*omega_sigma_2*omega_mean/T_omega)).*dw
        omegap[:,t+1] = omegap[:,t+1].*(omegap[:,t+1].>0) #enforcing positivity

        #stepping the decorrelation times
        t_decorr_p = t_decorr_p.-dt;
        t_decorr_m = t_decorr_m.-dt;
        #sets index of those to be renewed with 0 - which doesn't correspond to any particle
        phi_pm[1,(t_decorr_p.<=0)] .= 0 #the t_decorrs are 2d for some reason
        phi_pm[2,(t_decorr_m.<=0)] .= 0

        #split into cells, compute centres/targets, run ODE step
        eval_by_cell!(function (i,j,cell_particles)
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
            return nothing
        end, x_pos[:,t], y_pos[:,t], space_cells)

        phi_c = 0.5.*(phip[:,phi_pm[1,:],t]+phip[:,phi_pm[2,:],t])
        diffusion = zeros(2,np)
        diffusion[1,:] = (phip[1,:,t]-phi_c[1,:]).*(exp.(-c_phi.*0.5.*omegap[:,t].*dt).-1.0)
        diffusion[2,:] = (phip[2,:,t]-phi_c[2,:]).*(exp.(-c_phi.*0.5.*omegap[:,t].*dt).-1.0)
        reaction = zeros(2,np) # body reaction
        # reaction = dt.*(reaction).*exp.(c_phi.*0.5.*omegap[:,t].*dt) #integration of reation term to match diffusion scheme, uncomment if reaction !=0
        dphi = (diffusion .+ reaction)

        # ensuring mean 0 change
        # generating a random orthonormal basis
        # is 2-d so genrate a random unit vector from an angle and proceed based
        # on that
        angle = 2*pi*rand(T, 1)[1];
        e_1 = [cos(angle),sin(angle)]
        handedness = sb.sample([-1,1],1)[1] #randomly choose betwen left or right handed system
        e_2 = handedness*[e_1[2],-e_1[1]]
        T_mat=zeros(2,2)
        T_mat[:,1] = e_1  #coord transform matrix
        T_mat[:,2] = e_2
        dphi = T_mat\dphi  #transform to new coords
        #performing adjustment to mean 0
        corr_factor = zeros(T, 2,np)
        
        eval_by_cell!(function (i,j,cell_particles)
            for phi_i=1:2
                phi_mean = mean(dphi[phi_i,cell_particles])
                if phi_mean != 0 #isn't true for empty cells
                    cell_points_pos = dphi[phi_i,cell_particles].>0
                    cell_points_neg = dphi[phi_i,cell_particles].<0
                    phi_pos_mean = mean(cell_points_pos.*dphi[phi_i,cell_particles])
                    phi_neg_mean = mean((cell_points_neg).*dphi[phi_i,cell_particles])
                    if phi_mean>0
                        corr_factor[phi_i,cell_particles]=.-cell_points_pos*(phi_neg_mean./phi_pos_mean) + (1 .- cell_points_pos)
                    else
                        corr_factor[phi_i,cell_particles]=.-(cell_points_neg)*(phi_pos_mean./phi_neg_mean) + (1 .- cell_points_neg)
                    end
                end
            end
            return nothing
        end, x_pos[:,t], y_pos[:,t], space_cells)

        dphi = corr_factor.*dphi
        dphi = T_mat*dphi #return to old coords
        phip[:,:,t+1] = phip[:,:,t]+dphi
        if !(initial_condition == "triple delta")
            phip[:,:,t+1] = phip[:,:,t+1].*(phip[:,:,t+1].>0) #forcing positive concentration
        end
        #bcs give modification for phi based on interaction between t=t and t=t+1
        set_phi_as_ic_up1!(phip,t+1, bc_interact[:,t,3] .| bc_interact[:,t,4])#set incoming particles to same as ic
        phip_diff = phip[:,:,t+1]
        bc_absorbtion!(phip, bc_interact[:,t,2], bc_params, t+1, precomp_P) #currently only reacting on bottom bc
        phip_diff -= phip[:,:,t+1]
        try
            eval_by_cell!((i,j,cell_particles)-> (if i==1;assign_f_phi_cell!(y_0_flux,phip_diff[:,cell_particles[bc_interact[cell_particles,t,2]]],psi_mesh,i,j,t)
            ; end;return nothing) , x_pos[:,t+1], y_pos[:,t+1], space_cells)
        catch
            
            rethrow()
        end

        assign_f_phi!(f_phi,phip[:,:,t], x_pos[:,t], y_pos[:,t], psi_mesh, space_cells,t)
        # print(maximum(phip[:,:,t]),' ')
    end
    verbose && println("end")
    return nothing
end

function make_f_phi_no_PSP!(f_phi::Array{T,5},x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::Array{T,2}, bc_interact::Array{Bool,3}, initial_condition::String, psi_mesh::PsiGrid{T}, space_cells::CellGrid{T}, bc_params::BCParams{T}, verbose::Bool=false) where T<:AbstractFloat
    np, nt = size(x_pos)
    nt-=1
    
    phip = zeros(T, (2, np, nt+1)) #scalar concentration at these points

    if initial_condition == "Uniform phi_1"
        set_phi_as_ic_up1!(phip,1)
    elseif initial_condition == "triple delta"
        set_phi_as_ic_td!(phip,1)
    elseif initial_condition == "2 layers"
        set_phi_as_ic_2l!(phip,y_pos[:,1],space_cells,1)
    elseif initial_condition == "double delta"
        set_phi_as_ic_dd!(phip,1)
    else
        throw(ArgumentError("Not a valid intitial condition"))
    end
    assign_f_phi!(f_phi,phip[:,:,1], x_pos[:,1], y_pos[:,1], psi_mesh, space_cells,1)
    
    for t in 1:nt
        verbose && print(t,' ')
        bc_absorbtion!(phip,bc_interact[:,t,2],turb_k_e[bc_interact[:,t,2],t],bc_params,t) #currently only reacting on bottom bc
        phip[:,:,t+1] = phip[:,:,t]
        assign_f_phi!(f_phi,phip[:,:,t+1], x_pos[:,t+1], y_pos[:,t+1], psi_mesh, space_cells,t+1)
    end
    verbose && println("end")
    return nothing
end

#constant turb_k_e
function make_f_phi_no_PSP!(f_phi::Array{T,5},x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::T, bc_interact::Array{Bool,3}, initial_condition::String, psi_mesh::PsiGrid{T}, space_cells::CellGrid{T}, bc_params::BCParams{T}, verbose::Bool=false) where T<:AbstractFloat
    np, nt = size(x_pos)
    nt-=1
    precomp_P = min.(bc_params.bc_k.*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)

    phip = zeros(T, (2, np, nt+1)) #scalar concentration at these points

    if initial_condition == "Uniform phi_1"
        set_phi_as_ic_up1!(phip,1)
    elseif initial_condition == "triple delta"
        set_phi_as_ic_td!(phip,1)
    elseif initial_condition == "2 layers"
        set_phi_as_ic_2l!(phip,y_pos[:,1],space_cells,1)
    elseif initial_condition == "double delta"
        set_phi_as_ic_dd!(phip,1)
    else
        throw(ArgumentError("Not a valid intitial condition"))
    end
    assign_f_phi!(f_phi,phip[:,:,1], x_pos[:,1], y_pos[:,1], psi_mesh, space_cells,1)
    
    for t in 1:nt
        verbose && print(t,' ')
        bc_absorbtion!(phip, bc_interact[:,t,2], bc_params, t, precomp_P) #currently only reacting on bottom bc
        phip[:,:,t+1] = phip[:,:,t]
        assign_f_phi!(f_phi,phip[:,:,t+1], x_pos[:,t+1], y_pos[:,t+1], psi_mesh, space_cells,t+1)
    end
    verbose && println("end")
    return nothing
end

function make_f_phi_no_PSP_inflow_record_flux!(f_phi::Array{T,5},y_0_flux::Array{T,5},x_pos::Array{T,2},y_pos::Array{T,2}, turb_k_e::T, bc_interact::Array{Bool,3}, initial_condition::String, psi_mesh::PsiGrid{T}, space_cells::CellGrid{T}, bc_params::BCParams{T}, verbose::Bool=false) where T<:AbstractFloat
    np, nt = size(x_pos)
    nt-=1
    precomp_P = min.(bc_params.bc_k.*sqrt.(bc_params.B.*pi./(bc_params.C_0.*turb_k_e)),1)

    phip = zeros(T, (2, np, nt+1)) #scalar concentration at these points

    if initial_condition == "Uniform phi_1"
        set_phi_as_ic_up1!(phip,1)
    elseif initial_condition == "triple delta"
        set_phi_as_ic_td!(phip,1)
    elseif initial_condition == "2 layers"
        set_phi_as_ic_2l!(phip,y_pos[:,1],space_cells,1)
    elseif initial_condition == "double delta"
        set_phi_as_ic_dd!(phip,1)
    else
        throw(ArgumentError("Not a valid intitial condition"))
    end
    
    for t in 1:nt
        verbose && print(t,' ')
        set_phi_as_ic_up1!(phip,t, bc_interact[:,t,3] .| bc_interact[:,t,4])#set incoming particles to same as ic
        phip_diff = phip[:,:,t]
        bc_absorbtion!(phip, bc_interact[:,t,2], bc_params, t, precomp_P) #currently only reacting on bottom bc
        phip_diff -= phip[:,:,t]
        phip[:,:,t+1] = phip[:,:,t]
        assign_f_phi!(f_phi,phip[:,:,t], x_pos[:,t], y_pos[:,t], psi_mesh, space_cells,t)
        eval_by_cell!((i,j,cell_particles)-> (if i==1;assign_f_phi_cell!(y_0_flux,phip_diff[:,cell_particles],psi_mesh,i,j,t)
        ; end;return nothing) , x_pos[:,t], y_pos[:,t], space_cells)
    end
    verbose && println("end")
    return nothing
end
end