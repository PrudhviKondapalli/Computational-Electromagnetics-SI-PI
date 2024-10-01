%%%%%% Parameters Setup
mu=4*pi*1e-7; eps=8.853e-12; c0=1/sqrt(mu*eps);
f=10e9; lambda=c0/f; omega=2*pi*f;
dx = lambda/20;
dt = (dx/(c0));
Nx = 100; % # of FDTD cells
Cmu = dt/(mu*dx);
Cep = dt/(eps*dx);
Ez = zeros(Nx+1,1);
Hy = zeros(Nx,1);
%%%%%% Time marching = Leap-frog loop
for n=1:300 
    wt=omega*n*dt;
    %%Inorder to implement ABC at right boundary, saving Hy value at Nx
    oldHy = Hy(Nx);
    for i=1:Nx
        Hy(i) = Hy(i)+Cmu*(Ez(i+1)-Ez(i));
    end
    for i=2:Nx % i=1,Nx+1 omitted to simulate PEC boundaries
        Ez(i) = Ez(i)+Cep*(Hy(i)-Hy(i-1));
    end
    %%We have to update the Ez at the boundary based on the saved Hy(Nx)
    Ez(Nx+1) = Ez(Nx+1) + Cep*(oldHy - Hy(Nx));

    %%Hard source
    if wt < 2*pi
        Ez(Nx/2+1) = 1/32*(10-15*cos(wt)+6*cos(2*wt)-cos(3*wt)); % Hard Source
    else
        Ez(Nx/2+1)=0;
    end

%     if wt < 2*pi
%         Ez(Nx/2+1) = Ez(Nx/2+1) + 1/32*(10-15*cos(wt)+6*cos(2*wt)-cos(3*wt)); % Soft Source
%     end

subplot(211),plot(Ez(1:Nx)); axis([1 100 -1 1]);title('Ez')
subplot(212),plot(Hy(1:Nx)); axis([1 100 -6e-3 6e-3]); title('Hy'); pause(0.05)
end
% plot(Ez(1:Nx))

