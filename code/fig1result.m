% Author: Reza Ghasemi Alavicheh [r.ghasemi.reza@gmail.com]

clc;clear

%% random parameters
[X,C] = materni(cluster,"on");
pd = makedist('Nakagami','mu',1,'omega',1);
h = random(pd,1,2);

%% parameters
L =11;
cluster = 4 ;
iter = 1000;
R_B = 1700;
W = 100:100:500;
power = 1:1:5;
[W_m, power_m] = meshgrid(W,power);
pr = zeros(length(power),length(power_m));
for i = 1:1:length(power_m)
    for j = 1:1:length(W_m)
        [probOFasso,~] = func_s2_without_BWconstraint(iter,L,R_B,W(i),power(j),X,h,C,cluster);
        pr(j,i) = probOFasso;
    end
end
figure
[C,h] = contour(W_m,power_m,pr,LineWidth=1.6);
clabel(C,h);
xlabel('BandWidth(MHz)')
ylabel('Transmitted Power (w)')
grid on;
legend('success probability of association')