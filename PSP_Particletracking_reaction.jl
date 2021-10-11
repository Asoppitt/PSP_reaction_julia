using Random, Distributions
import StatsBase as sb
import LinearAlgebra as la
import Statistics as st

Random.seed!(12345) #setting a seed

n=9 #a convinent method of increascing resolution while maintaining 
    # resolution ratios and number of particles per cell

np = n*n*1000 # number of particles
dt = 0.01  # time step
nt = 20   # number of time steps
length_domain = 320 #length of periodic element
height_domain = 320
phi_domain = [0,1.2]
omega_shape_init = 0.25
omega_sigma_2 = 0.25
psi_partions_num = 20 #number of partitions in 1 state space direction
c_phi = 1.2
c_t = 2
u_max = 0.06
function u_mean(y)
    -4*u_max.*y.*(y.-height_domain)./height_domain^2
end
C_0 = 1.2  # rate parameter for velocity change
B=(1.5*C_0) #for reducing the langevin to a standard form - part of boundary conditions implementaion from Erban and Chapman 06
x_res=n #number of cells in x dim
y_res=ceil(Int, n/2)

bc_k=1
omega_mean=10
turb_e_init=1
T_omega = 1/(omega_mean); #approximation to frequency timescale

#ititialising position
xp = length_domain*rand(Float64, (np,nt+1)) # position of particle in x-direction
yp = height_domain*rand(Float64, (np,nt+1)) # position of paticle in y-direction
phip = zeros((2, np, nt+1)) #scalar concentration at these points

#ititialising at 1 for phi_1 for 1-D testing
phip[2,:,1] = 0.001*randn((np,1)) #pdf can't find zeros
phip[1,:,1] .= 1;

omegap = zeros(np,nt+1) #turbulence frequency
omega0_dist = Gamma(omega_mean/omega_shape_init,omega_shape_init)
omegap[:,1] = rand(omega0_dist, np)

x_edges = LinRange(0,length_domain,x_res+1)
y_edges = LinRange(0,height_domain,y_res+1)

#array for storing indicies of paired particles for psp
phi_pm = zeros(Int, 2, np)

function assign_pm(particles::Vector{Int}, cell_points::Vector{Int}, time_index::Int)
    # might be worth a strategy that removes tested pairs, can't see how to make it not require a biiiig temp variable though 
    for particle in particles
        count = 0
        while true
            count = count+1
            try_pm = sb.sample(cell_points,2,replace=false, ordered=true) #ordering added to reduce sample space as order for pm pairs is irrlevent
            test_dot = la.dot((phip[:,try_pm[1],time_index]-phip[:,particle,time_index]),(phip[:,try_pm[2],time_index]-phip[:,particle,time_index]))
            if  test_dot<= 0
                phi_pm[:,particle] = try_pm
                break
            end
        end
    end
end

function assign_pm_single(particles::Vector{Int}, cell_points::Vector{Int}, index_tbc::Int, time_index::Int)
    known_index = 1-(index_tbc-1)+1 #flips 1 to 2 and 2 to 1
    for particle in particles
        count = 0
        while true
            count = count+1
            try_pm = sb.sample(cell_points,1,replace=false) 
            test_dot = la.dot((phip[:,try_pm[1],time_index]-phip[:,particle,time_index]),(phip[:,known_index,time_index]-phip[:,particle,time_index]))
            if  test_dot<= 0
                phi_pm[index_tbc,particle] = try_pm[1]
                break
            end
        end
    end
end

psi_1 = LinRange(phi_domain[1], phi_domain[2], psi_partions_num+1) #defines boundaries for cells in psi space
psi_2 = LinRange(phi_domain[1], phi_domain[2], psi_partions_num+1)

f_phi = zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1) # histogram of each cell at a given time 

function assign_f_phi_cell(cell_points::Vector{Int}, cell_row::Int, cell_column::Int, t_index::Int)
    "use if cell_points has alreday been determined"
    for psi1_i=1:psi_partions_num
        in_1 = psi_1[psi1_i].<phip[1,cell_points,t_index].<psi_1[psi1_i+1]
        for psi2_i=1:psi_partions_num
            in_2 = psi_2[psi2_i].<phip[2,cell_points,t_index].<psi_2[psi2_i+1]
            f_phi[psi1_i, psi2_i, cell_row, cell_column, t_index] = sum(in_1.&in_2)
        end
    end
    f_phi[:, :, cell_row, cell_column, t_index] = f_phi[:, :, cell_row, cell_column, t_index]./length(cell_points)
end 

function assign_f_phi(t_index::Int)
    for i in 1:y_res
        in_y = y_edges[i].<yp[:,t_index].<y_edges[i+1]
        for j in 1:x_res
            in_x = x_edges[j].<xp[:,t_index].<x_edges[j+1]
            cell_particles = findall(in_x.&in_y)
            assign_f_phi_cell(cell_particles, i, j, t_index)
        end 
    end
end 

#initialising phi_pm and f_phi
for i in 1:y_res
    in_y = y_edges[i].<yp[:,1].<y_edges[i+1]
    for j in 1:x_res
        in_x = x_edges[j].<xp[:,1].<x_edges[j+1]
        cell_particles = findall(in_x.&in_y)
        assign_pm(cell_particles,cell_particles,1)
        assign_f_phi_cell(cell_particles, i, j, 1)
    end 
end

#time stamp until new p/m bound found, needed to ensure particles are
#decorreltaed
t_decorr_p = 1 ./(c_t.*omegap[phi_pm[1,:],1])
t_decorr_m = 1 ./(c_t.*omegap[phi_pm[2,:],1])

k=zeros(np)
k= k.+turb_e_init

#intitial vaules of velocity, maintainingconsitancy with energy
uxp = randn(np).*sqrt.(2/3*k)
uyp = randn(np).*sqrt.(2/3*k)

function bc_absorbtion(abs_points::Vector{Int}, n_abs::Int, t_index::Int)
    abs_k = bc_k.*phip[:, abs_points, t_index]
    #K for Erban and Chapman approximation 
    P = zeros(n_abs,2)
    P[:,1] = abs_k[1,:].*sqrt.(B.*pi./(C_0.*k[abs_points]))
    P[:,2] = abs_k[2,:].*sqrt.(B.*pi./(C_0.*k[abs_points]))
    xi = rand(n_abs,2)
    phip[1, abs_points[P[:,1].>xi[:,1]], t_index] .= 0.0001
    phip[2, abs_points[P[:,2].>xi[:,2]], t_index] .= 0.0001
    
end

for t in 1:nt
    #Using the SLM model for particle location
    for i in 1:y_res
        in_y = y_edges[i].<yp[:,t].<y_edges[i+1]
        for j in 1:x_res
            in_x = x_edges[j].<xp[:,t].<x_edges[j+1]
            cell_particles = findall(in_x.&in_y)
            k[cell_particles].=0.5*(st.mean(uxp[cell_particles].^2)+st.mean(uyp[cell_particles].^2))*1.5 #turb_e_init;
        end
    end
    global uxp = uxp+(-0.5*B*omega_mean*uxp)*dt+randn(np).*sqrt.(C_0.*k.*omega_mean.*dt); 
    global uyp = uyp+(-0.5*B*omega_mean*uyp)*dt+randn(np).*sqrt.(C_0.*k.*omega_mean.*dt); 
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

    # bc_absorbtion(mag, dim_mag, t) disabled to match paper

    # Reflection at lower boundary y<0
    mag = findall(yp[:,t+1].<=0) # index of particle with yp>height_domain
    dim_mag = size(mag) # dimension of array "mag"

    y_mag_succ = yp[mag,t+1] # yp at time t+1 corresponding to the index "mag"

    ypr_mag = - y_mag_succ  # yp at time t+1 of the reflected particle

    yp[mag,t+1]= ypr_mag #replacement of yp>1 with yp of reflected particle
    uyp[mag] = -uyp[mag] #reflecting velocity

    bc_absorbtion(mag, dim_mag[1], t)

    #bc at end (y=length_domain) of domain
    end_indicies = xp[:,t+1].>=length_domain #index of particle with xp>length
    dim_end = sum(end_indicies)

    end_x = xp[end_indicies,t+1]
    xpr_end = end_x .- length_domain #shifting particles back to begining
    xp[end_indicies,t+1] = xpr_end #replacing x coords

    #resetting concentration of "new" particles
    phip[2,end_indicies,t+1] = 0.001*randn(dim_end,1) #pdf can't find zeros
    phip[1,end_indicies,t+1] .= 1

    #resetting energy of "new" particles
    uxp[end_indicies] = randn(dim_end,1)*sqrt(2/3*turb_e_init)
    uyp[end_indicies] = randn(dim_end,1)*sqrt(2/3*turb_e_init)
    k[end_indicies] .= turb_e_init

    #bc at start (y=0) of domain
    start_indicies = xp[:,t+1].<=0 #index of particle with xp>length
    dim_start = sum(start_indicies)

    xpr_start = -xp[start_indicies,t+1] #shifting particles back to begining
    xp[start_indicies,t+1] = xpr_start #replacing x coords
    uxp[start_indicies] = - uxp[start_indicies]

    #resetting concentration of "new" particles
    phip[2,start_indicies,t+1] = 0.001*randn(dim_start,1) #pdf can't find zeros
    phip[1,start_indicies,t+1] .= 1

    #resetting energy of "new" particles
    uxp[start_indicies] = randn(dim_start,1)*sqrt(2/3*turb_e_init)
    uyp[start_indicies] = randn(dim_start,1)*sqrt(2/3*turb_e_init)
    k[start_indicies] .= turb_e_init


    #PSP begins here
    #E-M solver for omega
    dw = sqrt(dt).*randn(np) #random draws
    omegap[:,t+1] = omegap[:,t]-(omegap[:,t].-omega_mean)./T_omega.*dt + sqrt.(omegap[:,t].*(2*omega_sigma_2*omega_mean/T_omega)).*dw
    omegap[:,t+1] = omegap[:,t+1].*(omegap[:,t+1].>0) #enforcing positivity

    #stepping the decorrelation times
    global t_decorr_p = t_decorr_p.-dt;
    global t_decorr_m = t_decorr_m.-dt;
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
            #update pairs to ensure all are within the same bounds
            p_nin = .!(in.(phi_pm[1,cell_particles],Ref(cell_particles)))
            m_nin = .!(in.(phi_pm[2,cell_particles],Ref(cell_particles)))
            pm_nin = p_nin .& m_nin
            p_nin = xor.(p_nin , pm_nin)
            m_nin = xor.(m_nin , pm_nin)
            assign_pm(cell_particles[pm_nin], cell_particles,t)
            assign_pm_single(cell_particles[p_nin], cell_particles, 1, t)
            assign_pm_single(cell_particles[m_nin], cell_particles, 2, t)
            t_decorr_p[cell_particles[pm_nin.&p_nin]] = 1 ./(c_t.*omegap[phi_pm[1,cell_particles[pm_nin.&p_nin]],1])
            t_decorr_m[cell_particles[pm_nin.&p_nin]] = 1 ./(c_t.*omegap[phi_pm[2,cell_particles[pm_nin.&m_nin]],2])
        end 
    end
    phi_c = 0.5.*(phip[:,phi_pm[1,:],t]+phip[:,phi_pm[2,:],t])
    diffusion = zeros(2,np)
    diffusion[1,:] = c_phi.*0.5.*omegap[:,t].*(phip[1,:,t]-phi_c[1,:])
    diffusion[1,:] = c_phi.*0.5.*omegap[:,t].*(phip[2,:,t]-phi_c[2,:])
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
                cell_points_pos = dphi[phi_i,cell_particles].>0
                phi_mean = mean(dphi[phi_i,cell_particles])
                phi_pos_mean = mean(cell_points_pos.*dphi[phi_i,cell_particles])
                phi_neg_mean = mean((1 .- cell_points_pos).*dphi[phi_i,cell_particles])
                if phi_mean>0
                    corr_factor[phi_i,cell_particles]=.-cell_points_pos*(phi_neg_mean./phi_pos_mean) + (1 .- cell_points_pos)
                else
                    corr_factor[phi_i,cell_particles]=.-(1 .- cell_points_pos)*(phi_pos_mean./phi_neg_mean) + cell_points_pos
                end
            end
        end
    end
    dphi = corr_factor.*dphi;
    dphi = T*dphi #return to old coords

    phip[:,:,t+1] = phip[:,:,t]+dphi
    phip[:,:,t+1] = phip[:,:,t+1].*(phip[:,:,t+1].>0) #forcing positive concentration

    assign_f_phi(t)
    println(maximum(f_phi[:,:,:,:,t]))
    println(minimum(f_phi[:,:,:,:,t]),'\n')

end

write("new_reaction_aved_5_flux", f_phi)
print("Success")