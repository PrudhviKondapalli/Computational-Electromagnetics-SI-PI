%%%%%% Parameters Setup 
mu=4*pi*1e-7; eps=8.853e-12; c0=1/sqrt(mu*eps); f=10e9; lambda=c0/f; omega=2*pi*f; dx = lambda/20;
dt1 = (dx/c0); dt2 = 0.9*(dx/c0); % Two different dt values
Nx = 100; % # of FDTD cells
Cmu1 = dt1/(mu*dx); Cmu2 = dt2/(mu*dx);
Cep1 = dt1/(eps*dx); Cep2 = dt2/(eps*dx);
Ez1 = zeros(Nx+1,1); Hy1 = zeros(Nx,1); % Arrays for case 1
Ez2 = zeros(Nx+1,1); Hy2 = zeros(Nx,1); % Arrays for case 2

figure;
subplot(2,1,1);
p1 = plot(Ez1(1:Nx), 'b-', 'LineWidth', 2);
hold on;
p2 = plot(Ez2(1:Nx), 'r--', 'LineWidth', 2);
axis([1 100 -1.5 1.5]);
title('Ez');
legend('dt = (dx/c0) at dx = lambda/20', 'dt = 0.9*(dx/c0) at dx = lambda/20');

subplot(2,1,2);
p3 = plot(Hy1(1:Nx), 'b-', 'LineWidth', 2);
hold on;
p4 = plot(Hy2(1:Nx), 'r--', 'LineWidth', 2);
axis([1 100 -6e-3 6e-3]);
title('Hy');
legend('dt = (dx/c0) at dx = lambda/20', 'dt = 0.9*(dx/c0) at dx = lambda/20');

%%%%%% Time marching = Leap-frog loop
for n=1:150
    wt1=omega*n*dt1;
    wt2=omega*n*dt2;
    for i=1:Nx
        Hy1(i) = Hy1(i)+Cmu1*(Ez1(i+1)-Ez1(i));
        Hy2(i) = Hy2(i)+Cmu2*(Ez2(i+1)-Ez2(i));
    end
    for i=2:Nx % i=1,Nx+1 omitted to simulate PEC boundaries
        Ez1(i) = Ez1(i)+Cep1*(Hy1(i)-Hy1(i-1));
        Ez2(i) = Ez2(i)+Cep2*(Hy2(i)-Hy2(i-1));
    end
    if wt1< pi
        Ez1(Nx/2) = 1; %Hard source pulse
    else
        Ez1(Nx/2) = 0;
    end
    if wt2< pi
        Ez2(Nx/2) = 1; %Hard source pulse
    else
        Ez2(Nx/2) = 0;
    end
    
    % Update plots
    set(p1, 'YData', Ez1(1:Nx));
    set(p2, 'YData', Ez2(1:Nx));
    set(p3, 'YData', Hy1(1:Nx));
    set(p4, 'YData', Hy2(1:Nx));
    drawnow;
    pause(0.001); % Adjust the pause duration for desired animation speed
end