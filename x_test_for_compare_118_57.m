tau = 0.3;
% % % Constructing the two dataset given by the Ashley M Hou et.al 2020
for ii = 1:nt_r_x
A = randraw('uniform', [-1, 1], [length(const.loadbuses),length(const.loadbuses)]);
Abar = A*A';
for k = 1:length(const.loadbuses)
    for j = 1:length(const.loadbuses)
        sig(k,j,ii) = tau*(Abar(k,j)/sqrt(Abar(k,k)*Abar(j,j)))*sqrt(const.d_hat(k)*const.d_hat(j));
    end
end
xi(:,ii) = mvnrnd(zeros(length(const.loadbuses),1),sig(:,:,ii))';
end
x_test = xi';