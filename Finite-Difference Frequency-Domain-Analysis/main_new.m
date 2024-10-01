%% ASSIGNMENT - 6
% Author: PRUDHVI KONDAPALLI
% Email: Prudhvi.kondapalli@colorado.edu
% Credits : Prof. Mohammed Hadi
% Course: ECEN 5524 Principles of Computational Electromagnetics for SI&PI

%% Question - 1
% Write a code to complete the shielded microstrip problem
% to obtain the first TM01 mode. 
% height(h) = width/2 = 2mm
% Width(w) = 4mm
% Dielectric constant(Er) = 4
% a = 40h = 80mm
% b = 2.5h = 5mm
% Mh(Into how many nodes(cells) is the height of substrate divided into)
% Mh should be taken greater than 10 and should be even
% Assume Frequency is 2 GHz
% Width of PEC Box = a = 40h = 80mm
% Height of PEC Box = b = 2.5h = 5mm
% Assume Mh is taken as 20 cells, then Mw = 2Mh = 40 cells
% Ma = 20Mh = 400 cells
% Mb = 10Mh = 200 cells

clc;
clear all

% Get the number cells that Height h should be divided into
Mh = input('Enter value of Mh :');
f = input('Enter value of Frequency :');
Mw = 2*Mh;
h = 2e-3;
w = 2*h;
Er = 4;
Ma = 40 * Mh;                   % As given in the question
Mb = 2.5 * Mh;                  % As given in the question

Uo = 4*pi*10^(-7);            % henry per mm
Eo = 8.854 * 10^ (-12);          % farad per mm

% Define three matrices - Ez, Sij, Kj
% Ez will be MaxMb matrix
Ez = zeros(Ma+1,Mb+1);

% Sij will be a MaxMb matrix
Sij = zeros(Ma+1,Mb+1);

% alpha will be Mh matrix
alpha = zeros(Mb+1);

% So as a intial exitation under the trace in the substrate can be 
% We can consider the i as the row just below the metal for this
% approximation
% For now I am doing all the i,j in the substrate as intial approximation. 
for i = 2:Ma
    for j = 2:(Mh)
        Ez(i,j) = sin((pi*j)/Mh);
        if j == 20
            Ez(i,j) = 0;
        end
    end
end
    
% Plot the initial Ez field distribution
figure(1);
mesh(Ez);
title('Initial Ez Field Distribution');
xlabel('Length PEC:5mm(b)');
ylabel('Height PEC:80mm(a)');
zlabel('Ez Field');
colorbar;

% First find K0
w = 2*pi*f;
ko_squared = (w^2) * Uo * Eo;
ko = sqrt(ko_squared);
betao_sqd = (w^2) * Uo * Eo;

% Now using this Ko, we have to find new Ez
% First compute all Alphas
for j = 2:Mb
    if(j > (Mh+1))
        alpha(j) = 4 - (ko_squared*((h/Mh)^2));
    end
    if (j == (Mh+1))
        x = ko_squared + ((betao_sqd)*((Er-1)/2));
        alpha(j) = 4 - (x*((h/Mh)^2));
    end
    if (j < (Mh+1))
        x = ko_squared + ((betao_sqd)*((Er-1)));
        alpha(j) = 4 - (x*((h/Mh)^2));
    end
end

% So now we have initial values, we need to calculate Ez over the entire
% domain atleast 3 to 5 times.
for a = 1:5
    for i = 2:Ma
        for j = 2:Mb
            Sij(i,j) = Ez(i+1,j) + Ez(i-1,j) + Ez(i,j+1) + Ez(i,j-1);
        end
    end
    % Update the Ez matrix using the new Sij values and alpha values
    for i = 2:Ma
        for j = 2:Mb
            Ez(i,j) = Sij(i,j) / alpha(j);
        end
    end
    % Remove all nodes on the metal to and put Ez value as 0
    for i = (((Ma/2)+1)-((Mw/2)+1)):(((Ma/2)+1)+((Mw/2)+1))
             Ez(i,(Mh+1)) = 0;
    end
end

% Assume
ko_diff = 1;
tolerance = 0.0001;
counter = 0;

%%while ko_diff >= tolerance
while counter < 1001
    counter = counter+1;
    ko_prev = ko;
    % Now to find the updated Ko, we need to find I1 and I2(I2a and I2b)
    I1 = 0;
    for i= 2:Ma
        for j = 2:Mb
         I1 = I1 + (Ez(i,j)*(Sij(i,j)-4*Ez(i,j)));
        end
    end

    I2 = 0;
    I2a = 0;
    I2b = 0;

    for i=2:Ma
        for j=2:Mb
            I2a = I2a + (Ez(i,j)^2);
        end
    end

    I2a = ((h/Mh)^2)*I2a;

    for i=2:Ma
        I2b = I2b + (1/2 * betao_sqd)*(Er-1)*((h/Mh)^2)*(Ez(i,Mh)^2);
        for j=1:Mh
            I2b = I2b + (betao_sqd)*(Er-1)*((h/Mh)^2)*(Ez(i,j)^2);
        end
    end

    % Now ko can be recalculated
    ko_squared = -((I1 + I2b) / I2a);
    ko = sqrt(ko_squared);

    % Calculate the difference between the previous and current ko
    ko_diff = abs(ko - ko_prev);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First compute all Alphas
    for j = 2:Mb
        if(j > (Mh+1))
            alpha(j) = 4 - (ko_squared*((h/Mh)^2));
        end
        if (j == (Mh+1))
            x = ko_squared + ((betao_sqd)*((Er-1)/2));
            alpha(j) = 4 - (x*((h/Mh)^2));
        end
        if (j < (Mh+1))
            x = ko_squared + ((betao_sqd)*((Er-1)));
            alpha(j) = 4 - (x*((h/Mh)^2));
        end
    end
    
    % Now we have to use over relaxation to distribute Ez and find new Sij
    R = 3;
    for a = 1:R 
        for i = 2:Ma
            for j = 2:Mb
                Sij(i,j) = Ez(i+1,j)+Ez(i-1,j)+Ez(i,j+1)+Ez(i,j-1);
            end
        end

    % Now use Sij and alphaj to find new Ez
    for i = 2:Ma
        for j = 2:Mb
            Ez(i,j) = (Sij(i,j))/alpha(j);
        end    
    end

    % Remove all nodes on the metal to and put Ez value as 0
        for i = (((Ma/2)+1)-((Mw/2)+1)):(((Ma/2)+1)+((Mw/2)+1))
                Ez(i,Mh+1) = 0;
        end
    end
end
counter

% Plot the final Ez field distribution
figure(2);
mesh(Ez);
title('Final Ez Field Distribution');
xlabel('Length PEC:5mm(b)');
ylabel('Height PEC:80mm(a)');
zlabel('Ez Field');
colorbar;

% Plot the contour of the Ez field distribution
figure(3);
contourf(Ez);
title('Contour of Ez Field Distribution');
xlabel('Length PEC:5mm(b)');
ylabel('Height PEC:80mm(a)');
colorbar;

% Calculate the energy in different regions
for i = 2:Ma
    for j = Mh+2:Mb
    Energy_air = (1/2)*Eo*(Ez(i,j).^2); %Air
    end
end
Energy_sub = 0;
Energy_sub_1 = 0;
Energy_sub_2 = 0;

for i = 2:Ma
    for j = 2:Mh
    Energy_sub = Energy_sub + (1/2)*Eo*Er*(Ez(i,j).^2); %Substrate
    end
end


for i = 2:(Ma-Mw)/2
    Energy_sub_1 = Energy_sub_1+(1/2)*Eo*Er*(Ez(i,Mh+1).^2); %Around trace
end

for i = (Ma+Mw)/2:Ma
    Energy_sub_2 = Energy_sub_2+(1/2)*Eo*Er*(Ez(i,Mh+1).^2);
end

Total_Energy = Energy_air + Energy_sub+Energy_sub_1 + Energy_sub_2;
Energy_sub
Energy_sub_1
Energy_sub_2
Energy_Percentage_Substrate = (Energy_sub/Total_Energy)*100

