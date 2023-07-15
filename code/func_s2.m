% Author: Reza Ghasemi Alavicheh [r.ghasemi.reza@gmail.com]

function [probOFasso,average_num,avg_sum_rate,avg_SCBSs_number,avg_bw_consumption] = func_s2(iter,L,R_B,W,power,X,h,C,cluster)

number_of_association = zeros(1,iter);
sing_prob = zeros(1,iter);
sum_rate_data = zeros(1,iter);
SCBSs_numbers_data = zeros(1,iter);
bw_consumption_data = zeros(1,iter);
for k=1:1:iter

%% SCBSs and UAVs positions with 4 cluster from kmenas
%% constants
omega_t = 10*log10(power); % in dB
f_c = 2e9; % 2GHz
alpha = 9.61; gamma = 2;
beta = 0.16;  c = 3e8;
mu_LOS = 1 ; % in dB
mu_NLOS = 20; % in dB
sigma_n2 = -125; % in dB
GAMMA_min = -5;  %in dB
r_SCBS = [20,40,60,80,100] + rand(1,5); % requested rate in Mbps
UAVs_number = cluster; 
omega_r_matrix = zeros(length(X),UAVs_number);
height = 0.3;

%% received power in all SCBSs from each UAV
LAMBDA_temp =zeros(length(X),UAVs_number);
s_temp = zeros(length(X),UAVs_number);
d_ij_temp = zeros(length(X),UAVs_number);
for j = 1:1:UAVs_number
    for i =1:1:length(X)
        d_ij = sqrt( (X(i,1)-C(j,1))^2 + (X(i,2)-C(j,2))^2);
        d_ij_temp(i,j) = d_ij;
        s = sqrt(d_ij^2 + height^2);
        s_temp(i,j) = s;
        k_0 =10*log10((4*pi*s*f_c/c)^gamma);
        phi =real(atan(height/d_ij)); 
        P_LOS = 1/(1+alpha*exp(-beta*(phi-alpha)));
        P_NLOS = 1-P_LOS;
        F = P_LOS*h(1) + P_NLOS*h(2);
        LAMBDA = k_0 + P_LOS*mu_LOS + P_NLOS*mu_NLOS - F; %in dB
        LAMBDA_temp(i,j) = LAMBDA;
        omega_r_matrix(i,j) = omega_t - LAMBDA; % in dB
    end
end
%% SINR in each SCBS from each UAV
GAMMA = zeros(length(X),UAVs_number);
omega_r_matrix_w = 10.^(omega_r_matrix./10);
sigma_n2_w = 10.^(sigma_n2./10);
for i = 1:1:length(X)
    for j=1:1:UAVs_number
        sum_omega = 0;
        for v = 1:1:UAVs_number %compute sum in SINR
            if v == j
                continue
            else
                sum_omega_temp = omega_r_matrix_w(i,v);
                sum_omega = sum_omega + sum_omega_temp;
            end
        end
        GAMMA(i,j) = real(10*log10(omega_r_matrix_w(i,j))-10*log10(sum_omega + sigma_n2_w));
    end
end
GAMMA_nondb = 10.^(GAMMA./10);
%% BW demanding of each SCBS
w = zeros(length(X),UAVs_number);
for j = 1:1:UAVs_number
    for i =1:1:length(X)
        w(i,j) = real(r_SCBS(randi(length(r_SCBS),1))/log2(1+GAMMA_nondb(i,j)));
    end
end
%% rate demanding of each SCBS
r = zeros(length(X),UAVs_number);
for j = 1:1:UAVs_number
    for i =1:1:length(X)
        r(i,j) = real(r_SCBS(randi(length(r_SCBS),1)));
    end
end

%% spectral Efficiency
% SE = r./w;
% SE_sort = sort(SE,1);

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%  Matrix A (Algorithm 1) %%%%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_1 = double(GAMMA > GAMMA_min); % compare SINR with minSINR and if TRUE put B_ij =1
A = B_1;
for i = 1:1:length(X)
    if sum(A(i,:)) > 1
        A(i,:) = 0;
        [~,J] = min(r(i,:)+w(i,:));
        A(i,J) = 1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UAV step
% SE_A = A.*SE ;
w_A = w.*A;
r_A = r.*A;
r_w_A = r_A+w_A; 
r_w_A_sort = sort(r_w_A,"ascend");
for j = 1:1:UAVs_number
    T_L = 0 ;
    T_B = 0;
    for kk =1:1:length(A)
        index_of_r_w_A = find(r_w_A(:,j) == r_w_A_sort(kk,j));
        if isempty(index_of_r_w_A)            
            continue
        end
        if r_w_A_sort(kk,j) == 0            
            continue
        end
        if (T_L < L) && (T_B < W)
            if T_B + w(index_of_r_w_A,j) <= W
                T_L = T_L + 1;
                T_B = T_B + w_A(index_of_r_w_A,j);
            else
                A(index_of_r_w_A,j) = 0;                
            end
        else
            A(index_of_r_w_A,j) = 0;
        end
    end
end

r_A = r.*A;
T_r =sum(sum(A.*r));
% T_r_0 = T_r;

% J_0 = 0;
flag = 0;
l =0;
sum_A_index = sum(A);
while sum(sum(A.*r)) > R_B
    [~,J] = max(sum_A_index);
    sum_A_index(J) =0;
    if sum(sum_A_index) == 0
        sum_A_index = sum(A);
        flag = 1;
        l = l +1;
    end
%     J_0 = J;
    for i = 1:1:length(X)
        r_A_J = r_A(:,J);
        if flag == 0
            if (r_A(i,J) == min((r_A_J(r_A_J>0))))
                A(i,J) = 0;
                T_r = T_r - r_A(i,J);
                T_L = T_L + 1;
                T_B = T_B + w_A(i,J);
            end
        else
            if (r_A(i,J) == min((r_A_J(r_A_J>r_SCBS(l)))))
                A(i,J) = 0;
                T_r = T_r - r_A(i,J);
                T_L = T_L + 1;
                T_B = T_B + w_A(i,J);
            end

        end
        
    end
end

prob =sum(sum(A))/length(X);
disp(prob);

number_of_association(1,k) = sum(sum(A));
sing_prob(1,k) = sum(sum(A))/length(A);
sum_rate_data(1,k) = sum(sum(A.*r));
SCBSs_numbers_data(1,k) = sum(sum(A));
bw_consumption_data(1,k) = sum(sum(A.*w))/cluster;
disp(k);

end
clc;
disp(length(X))
average_num = sum(number_of_association)/length(number_of_association);
disp(average_num)
probOFasso = sum(sing_prob)/length(number_of_association);
disp(probOFasso);
avg_sum_rate = sum(sum_rate_data/length(sum_rate_data));
avg_SCBSs_number = sum(SCBSs_numbers_data/length(SCBSs_numbers_data));
avg_bw_consumption = sum(bw_consumption_data/length(bw_consumption_data));
end