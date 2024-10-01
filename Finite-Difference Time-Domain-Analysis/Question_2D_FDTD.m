%%%%%% Parameters Setup
mu=4*pi*1e-7; eps=8.853e-12; c0=1/sqrt(mu*eps);
f=1e9; lambda=c0/f; omega=2*pi*f;

dx = lambda/40;
dy = lambda/40;

%Based on numerical stability criterion
%dt = 1*((sqrt(mu*eps))/((sqrt((2/dx)^2))));
dt = dx/(2*c0)/sqrt(2);

Nx = 100; % # of FDTD cells
Ny = 100;

Cmux = dt/(mu*dx);
Cmuy = dt/(mu*dy);
Cepx = dt/(eps*dx);
Cepy = dt/(eps*dy);

Ex = zeros(Nx,Ny+1);
Ey = zeros(Nx+1,Ny);

Hz = zeros(Nx,Ny);
%%%%%% Time marching = Leap-frog loop
for n=1:75 
    wt=omega*n*dt;
    for i=1:Nx
        for j=1:Ny
            Hz(i,j) = Hz(i,j)+(Cmuy*(Ex(i,j+1) - Ex(i,j)))-(Cmux*(Ey(i+1,j)-Ey(i,j)));
        end
    end
    for i=1:Nx
        for j = 2:Ny
            Ex(i,j) = Ex(i,j)+Cepy*(Hz(i,j)-Hz(i,j-1));
        end
    end

    for i=2:Nx
        for j = 1:Ny
            Ey(i,j) = Ey(i,j)-Cepx*(Hz(i,j)-Hz(i-1,j));
        end
    end




    if wt < 2*pi
        Ex((Nx/2)+1,(Ny/2)+1) = 1/32*(10-15*cos(wt)+6*cos(2*wt)-cos(3*wt)); % Hard Source
        %%Ey(Nx/2+1) = 1/32*(10-15*cos(wt)+6*cos(2*wt)-cos(3*wt)); % Hard Source
    else
        Ex(Nx/2+1,Ny/2+1)=0;
        %%Ey(Nx/2+1)=0;
    end

%     if wt < 2*pi
%         Ez(Nx/2+1) = Ez(Nx/2+1) + 1/32*(10-15*cos(wt)+6*cos(2*wt)-cos(3*wt)); % Soft Source
%     end

%subplot(311), plot(Ex(1:Nx)); axis([1 100 -100 100]); title('Ex');
%subplot(312), plot(Ey(1:Ny)); axis([1 100 -100 100]); title('Ey');
%subplot(313), plot(Hz(:)); axis([1 100 -6e-3 6e-3]); title('Hz'); pause(0.05);
mesh(Ex); pause(0.00005);

end
% plot(Ez(1:Nx))
