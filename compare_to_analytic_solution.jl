using Plots

n=9
nt = 20
psi_partions_num = 20
x_res=n #number of cells in x dim
y_res=ceil(Int, n/2)
in_data=zeros(psi_partions_num, psi_partions_num, y_res, x_res, nt+1)
read!("new_reaction_aved_5_flux",in_data)


psi_spacing = LinRange(0,1.2,psi_partions_num+1) #defines boundaries for cells in psi space
psi_spacing = 0.5 .* (psi_1[1:psi_partions_num]+psi_1[2:psi_partions_num+1])

phi_1_means = [sum([in_data[:,psi_2_index,i,j,t].*psi_spacing for psi_2_index in 1:psi_partions_num]) for i in 1:y_res, j in 1:x_res, t in 1:nt]
print(size(phi_1_means))