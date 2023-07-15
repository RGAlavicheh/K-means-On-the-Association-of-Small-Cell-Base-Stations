% Author: Reza Ghasemi Alavicheh [r.ghasemi.reza@gmail.com]

clc;clear
%% parameters
iter =2000;
L =1:1:15;
R_B = 1400;
W = 200;
power = 1.3;
cluster = 4;
[X,~] = materni(cluster,"off");
if length(X) >29 || length(X)<27
    clc;
    clear;
    figure_optional
    return
end
pd = makedist('Nakagami','mu',1,'omega',1);
h = random(pd,1,2);
for conventional = 1:1:2
    %%
    data_a = zeros(length(L),1);
    data_b = zeros(length(L),1);
    data_c = zeros(length(L),1);
    data_d = zeros(length(L),1);
    if conventional == 1
        [~,C] = materni(cluster,"off"); %kmeans in matern function
        LW = 1.6;
        C_kmeans = C;
    else

        [s_min] = s_computer(X);
        [C] = materni_s(s_min);
        C = C(1:3,:);
        cluster =length(C);
        LW =1;
        C_conventional = C;
    end

%
LS = ["-" "--"];
    if conventional == 1
        figure
    else
        color_v = ['b','b','c'];
        marker_v = ['.','.','.'];

    end
    if conventional == 2
        color_p = color_v(1);
        marker_p = marker_v(1);
    else
        color_p = 'r';
        marker_p = 'x' ;
    end

    for i = 1:1:length(L)
       [probOFasso,average_num,avg_sum_rate,~,avg_bw_consumption] = test_func_s2(iter,L(i),R_B,W,power,X,h,C,cluster);
        data_a(i,1) = average_num;
        data_b(i,1) = probOFasso;
        data_c(i,1) = avg_sum_rate;
        data_d(i,1) = avg_bw_consumption;

    end
    if conventional == 2
        color_p = color_v(2);
        marker_p = marker_v(2);

    else
        color_p = 'k';
        marker_p = 'v' ;
    end
    plot(L,data_b,Color=color_p,marker=marker_p ,MarkerEdgeColor='k',LineWidth=LW,LineStyle=LS(conventional));
     ylim([0 1]);
    xlim([1 15]);
    grid on;
    hold on;






legend("S2 k-means","S2 conventional");
    xlabel('L (number of links)','FontSize',12)
    ylabel("Avg.probability of association")
    title ("association comparison")

end
