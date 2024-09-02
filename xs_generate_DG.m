function [xs,xs_pen] = xs_generate_DG(data,rbus,nt_r,D,rated_cap,dist_type)
% 
% nt_r = 10000;
% dist_type = 'mix';
% data = case33bw;
for i = 1:D % D/2 becoz generating P and Q separately 
    switch dist_type
        case 'beta'
            xs_pen(:,i) = randraw('beta', [rated_cap*0.3 rated_cap 2 8], [1 nt_r])';
%             xs_pen((nt_r/2)+1:nt_r,i) = randraw('beta', [rated_cap*0.2 rated_cap 1 5], [1 nt_r/2])';
        case 'weibull'
            xs_pen(:,i) = randraw('weibull', [0 100/rated_cap rated_cap/1.5], [1 nt_r])';
        case 'laplace'
            xs_pen(:,i) = randraw('laplace',[ rated_cap*0.5 rated_cap*0.05], [1 nt_r])';
         case 'normal'
            xs_pen(:,i) = randraw('normal', [rated_cap/2 rated_cap/6], [1 nt_r])';
        case 'mix'
            xs_pen(:,i) = randraw('beta', [rated_cap*0.2 rated_cap 22 3], [1 nt_r])';
            % ----------------------
            xs_risk2(:,i) = randraw('laplace',[ rated_cap*0.8 rated_cap*0.15], [1 nt_r])';
            idx_sample1 = randi(nt_r,nt_r/2,1);
            xs_pen(idx_sample1,i) = xs_risk2(idx_sample1,i); %xs(idx_sample1,i+D/2) = xs_risk2(idx_sample1,i+D/2);
        case 'mix2'
            xs_pen(:,i) = randraw('beta', [rated_cap*0.2 rated_cap 22 3], [1 nt_r])';
            % ----------------------
            xs_risk2(:,i) = randraw('weibull', [0 100/rated_cap rated_cap/1.5], [1 nt_r])';
            idx_sample1 = randi(nt_r,nt_r/2,1);
            xs_pen(idx_sample1,i) = xs_risk2(idx_sample1,i); %xs(idx_sample1,i+D/2) = xs_risk2(idx_sample1,i+D/2);
      %--------------------------------------------------------------------------------------------- 
    end
xs = repmat(data.bus(rbus,3)',[nt_r,1]);
% xs = xs-xs_pen;

end
