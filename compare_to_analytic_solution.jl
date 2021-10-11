using Plots

n=9
nt = 20
psi_partions_num = 20
x_res=n #number of cells in x dim
y_res=ceil(Int, n/2)
in_data=zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
length_domain = 320 #length of periodic element
height_domain = 320
read!("new_reaction_aved_5_flux",in_data)


psi_spacing = LinRange(0,1.2,psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_spacing[1:psi_partions_num]+psi_spacing[2:psi_partions_num+1])

phi_1_means = [1/x_res.*sum([1/y_res.*sum(in_data[:,psi_2_index,i,j,t].*psi_spacing ) for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt]

y_integral_means = sum(phi_1_means, dims=1)[1,:,:]

plot(LinRange(0,height_domain,x_res),y_integral_means[:,end])