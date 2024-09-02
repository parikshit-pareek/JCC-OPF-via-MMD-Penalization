% Code for generating the normal coorelated distribution 
% Data Taken from 
%       K. Jhala, et. al, “The dominant influencer of volt-age fluctuation (divf) for power distribution system" 
%       IEEE Trans. Power Systems, vol. 34, no. 6, pp. 4847–4856, 2019.
% -----------------------------------------------------------------------------------------------------------

function [xs,xs_cell] = data_gen_divf(data,rbus,nt,N,D)
N_nt = nt+N;
% Three different covariance matrics from the paper
sig_1 = [1 -0.0447;-0.0447 0.2];
sig_2 = [2 -0.0894;-0.0894 0.4];
sig_3 = [3 -0.1342;-0.1342 0.6];

xp = repmat(data.bus(rbus,3)',[N_nt,1]);
xq = repmat(data.bus(rbus,4)',[N_nt,1]);

z_mu = zeros(2,1); % zero mean of 2-D
idx_1  = 1:3:34;
idx_2  = 2:3:34;
idx_3 = 3:3:34;
xs = [xp xq];

for i =1:length(idx_1)
    xs_1(:,:,i) = mvnrnd(z_mu,sig_1,N_nt)/100;
    xs(:,idx_1(i)) = xs(:,idx_1(i))-xs_1(:,1,i); % Placing real power
    xs(:,D/2+idx_1(i)) = xs(:,D/2+idx_1(i))-xs_1(:,2,i); % Placing reactive power
end

for i =1:length(idx_2)
    xs_2(:,:,i) = mvnrnd(z_mu,sig_2,N_nt)/100;
    xs(:,idx_2(i)) = xs(:,idx_2(i))-xs_2(:,1,i); % Placing real power
    xs(:,D/2+idx_2(i)) = xs(:,D/2+idx_2(i))-xs_2(:,2,i); % Placing reactive power
end


for i =1:length(idx_3)
    xs_3(:,:,i) = mvnrnd(z_mu,sig_3,N_nt)/100;
    xs(:,idx_3(i)) = xs(:,idx_3(i))-xs_2(:,1,i); % Placing real power
    xs(:,D/2+idx_3(i)) = xs(:,D/2+idx_3(i))-xs_3(:,2,i); % Placing reactive power
end

xs_cell = {xs_1, xs_2,xs_3};

end