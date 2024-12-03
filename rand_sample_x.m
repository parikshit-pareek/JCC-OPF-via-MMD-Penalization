function initx = rand_sample_x(N, D, xlimit)
% randomly sample N points in the space specified by xlimit
initx = zeros(N,D);

parfor j=1:N
initx(j,:) = rand(1,D).*(xlimit(2,:)-xlimit(1,:)) + xlimit(1,:);
end
% parfor j=N-50+1:N
% initx(j,:) = rand(1,D).*(xlimit(2,:)-xlimit(2,:)/1.1) + xlimit(2,:)/1.1;
% % end

%% Code for generating random samples sapertaly for P and Q

% for j = 1:D/2
%     initx(:,j) = rand(N,1)*(xlimit(2,j)-xlimit(1,j)) + xlimit(1,j);
% end
%  for j = D/2+1:D
%     initx(:,j) = rand(N,1)*(xlimit(2,j)-xlimit(1,j)) + xlimit(1,j);
end
% initx=sort(initx,1);