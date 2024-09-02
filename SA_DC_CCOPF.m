%% Large Sample Based DCOPF and Testing of Violation probability 

% clear
% clc

% format short g
% sys = '57bus';
% data = ext2int(pglib_opf_case57_ieee);
% const = ex_extract_ccDCOPF(data); %% Extract information from the MPC structure


SA_DCOPF.target_size = round((2/0.05)*(log(1/10^(-4))+2*const.ngen),0); % Number of samples used for SAA
tic
% solver_name = 'mosek';


[SA_DCOPF] = desired_dist_DCOPF_perunit_2(data,SA_DCOPF.target_size,sys,solver_name);


toc
%% Testing via Samples 
nt_r_x = 10000;
x_test = SA_DCOPF.xs_pen; 
g_x = SA_DCOPF.g;
parfor j= 1:nt_r_x
[A1_test(:,:,j),Ao_test(:,j)] = dcopf_a1ao(x_test(j,:)',data);
 F_x(:,j)   = A1_test(:,:,j)*g_x + Ao_test(:,j);
end
 
F_x > 10^(-3);
SA_DCOPF.vio_individual = sum(ans,2);
SA_DCOPF.eps_individual = SA_DCOPF.vio_individual/nt_r_x;
[SA_DCOPF.max_eps_individual,SA_DCOPF.max_c] = max(SA_DCOPF.eps_individual);


SA_DCOPF.vio_joint = sum(ans);
SA_DCOPF.vio_idx = SA_DCOPF.vio_joint >=1;
SA_DCOPF.eps_joint = (sum(SA_DCOPF.vio_idx)/nt_r_x);




