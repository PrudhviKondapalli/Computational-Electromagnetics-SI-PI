
clear all
clc
%Increment in the steps of 0.0001 from 0 to 1
x = 0 : 0.2 : 1;
%Assign all the expressions to the functions
N1 = [0,0.4673,0.90105,1.30105,1.6673,2];
N2= [0,-0.2378,0.8240,2.2560,2.9920,2];
N3=[0,0.4208,0.9336,1.5924,2.3928,2];
%Plot the experssions and annotate them properly
plot(x,N1, '-',LineWidth=0.5);
hold on
plot(x, N2, '-.',LineWidth=2);
hold on
plot(x, N3, '--', LineWidth=2);
hold on
xline(1, '.', Linewidth=2)
xlim([0 1])
ylim([-1 5])
ylabel("Phi value at intervals")
xlabel("X values")
title("Boundary problem solution for poissons equation")
legend('At N=1','At N=2', ...
    'At N=3');