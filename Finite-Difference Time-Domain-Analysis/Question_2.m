%%%%%% Parameters Setup
mu = 4*pi*1e-7;
eps = 8.853e-12;
c0 = 1/sqrt(mu*eps);

f = 1e9;
lambda = c0/f;
omega = 2*pi*f;

dx = lambda/40;
dy = lambda/40;

dt = 0.9/(c0*sqrt(1/(dx^2) + 1/(dy^2)));

Nx = 100; % # of FDTD cells in x-direction
Ny = 100; % # of FDTD cells in y-direction

Cmu_x = dt/(mu*dx);
Cmu_y = dt/(mu*dy);
Cep_x = dt/(eps*dx);
Cep_y = dt/(eps*dy);

Ez = zeros(Nx+1, Ny+1);
Ex = zeros(Nx+2, Ny+1); % Extending Ex by one cell in x-direction
Ey = zeros(Nx+1, Ny+2); % Extending Ey by one cell in y-direction
Hz = zeros(Nx+1, Ny+1);

%%%%%% Time marching = Leap-frog loop
for n = 1:300
    wt = omega*n*dt;

    % Update Ex
    for i = 2:Nx+1
        for j = 2:Ny
            Ex(i, j) = Ex(i, j) + Cep_y*(Hz(i, j) - Hz(i, j-1));
        end
    end

    % Update Ey
    for i = 2:Nx
        for j = 2:Ny+1
            Ey(i, j) = Ey(i, j) - Cep_x*(Hz(i, j) - Hz(i-1, j));
        end
    end

    % Update Hz
    for i = 2:Nx
        for j = 2:Ny
            Hz(i, j) = Hz(i, j) - Cmu_x*(Ey(i+1, j) - Ey(i, j)) + Cmu_y*(Ex(i, j+1) - Ex(i, j));
        end
    end

    % Hard source
    if wt < 2*pi
        Ez(Nx/2+1, Ny/2+1) = 1/32*(10 - 15*cos(wt) + 6*cos(2*wt) - cos(3*wt));
    else
        Ez(Nx/2+1, Ny/2+1) = 0;
    end

    % Visualization (optional)
    % Visualization (optional)
    subplot(221), plot(Ez(Nx/2+1, :)), axis([1 Ny+1 -1 1]), title('Ez')
    subplot(222), plot(Hz(Nx/2+1, :)), axis([1 Ny+1 -6e-3 6e-3]), title('Hz')
    subplot(223), plot(1:Nx+2, Ex(:, Ny/2+1)), axis([1 Nx+2 -1 1]), title('Ex')
    subplot(224), plot(1:Ny+2, Ey(Nx/2+1, :)), axis([1 Ny+2 -1 1]), title('Ey')
    pause(0.05)

end