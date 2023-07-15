function [s_min] = s_computer(X)
%% constants
f_c = 2e9; % 2GHz
alpha = 9.61; gamma = 2;
beta = 0.16;  c = 3e8;
mu_LOS = 1 ; % in dB
mu_NLOS = 20; % in dB
height = 0.3;
PL_max_db = 110;
[ppp] = conve();
ss = [];
for i = 1:1:length(X)
    for j = length(ppp)
        d_ij = sqrt( (X(i,1)-ppp(j,1))^2 + (X(i,2)-ppp(j,2))^2);
        %% compute d_min
        phi =real(atan(height/d_ij));
        P_LOS = 1/(1+alpha*exp(-beta*(phi-alpha)));
        P_NLOS = 1-P_LOS;
        t1 = sqrt(10.^((PL_max_db-P_LOS*mu_LOS - P_NLOS*mu_NLOS)./10));
        s = t1*(c/(4*pi*f_c));

        ss = [ss s];
    end
end

s_min = min(ss);
end