clc;
clf;
clear all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input parameters:

% Plate Initialization:
    x_length = .05;         % [m]
    y_length = .05;         % [m]
    nx = 100;                % number of grid points (Should be at least 1 per mm)
	ny = 100;                % number of grid points
% Plate Properties:
    T0 = 300; 			    % Initial ambient temperature [K]
	q = 840; 			    % laser power [W]
	k = 35; 			    % conductivity coefficient [W/m/K]
	cp = 800; 			    % specific heat [J/kg/K]
	rho = 7600.0; 		    % density [kg/m^3]
% Laser Initialization:
    HSP =[-0.02,0];         % Initial position of heat source from plate center (x,y) [m]
    vx = 20/1000; 	        % laser speed x [m/s] (from [mm/sec])
    vy = 0*50/1000; 	    % laser speed y [m/s]
% Simulation Options: 
    Max_Time = 2;           % Time simulated [s]
    Data_Points = 200;      % amount of data sets collected
    Temp_sensitivity = 50;  % pecent of maximum temperature shown on figure (sensitivity)
    Scale = 1;              % figure window scale
    Window_Pos = [0,0];     % Position of window from plate center (x,y) [m]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Error Check:
    if abs(HSP(1)) >= (x_length/2) || abs(HSP(2)) >= (y_length/2)
        disp('Error: Heat Source is not on the plate.')
        return
    end

[x,y] = meshgrid(linspace(-x_length/2,x_length/2,nx),linspace(-y_length/2,y_length/2,ny)); % [m]

Time = linspace(0,Max_Time,Data_Points);    % Defines all times data is taken
dt = Max_Time/Data_Points;                  % Defines the change in time for each calculation
T_pre=[]; T_diff=[]; Trial_Temps=[];        % Sets containers for Data
T= T0*ones(ny,nx);                          % Creates initial temperature matrix

for n = 1:1:length(Time) % Loop that creates all sets of data required (Dependent on time))
    T_pre = T;
    for i=2:ny-1     % Loop that Diffuses the heat through the plate (Dependent on time and position)
        for j=2:nx-1 % Gauss-Seidel method was used to produce diffusion equation
            T_diff(i,j)=(T_pre(i,j)+k*dt*(T(i-1,j)+T(i+1,j)+T(i,j-1)+T(i,j+1)))/(1+4*dt*k);
        end
    end
    for i=2:ny-1    % Loop that introduces the heat source (Dependent on time and position)
        for j=2:nx-1
            A = (rho*cp)^(1/2);
            B = 8*((pi*k*dt)^(3/2));
            C = (rho*cp)/(4*k*dt);
            D = ((x(i,j)-HSP(1))^2)+((y(i,j)-HSP(2)))^2;
            Q = q*dt;  %Laser strength [J]
            T(i,j) = T_diff(i,j) + ((Q*A)/B)*exp(-C*D); % Equation details seen in code notes
        end
    end
    Trial_Temps = [Trial_Temps,T]; %Collects all data in one array

    HSP(1) = HSP(1)+vx*dt;         %Updates the position of the heat source
    HSP(2) = HSP(2)+vy*dt;
end

Max_Temp = max(max(Trial_Temps));  %Finds the Maximum Temperature produced

figure(1)       % Creates a figure with annotations for the contours
annotation('textbox', [0.13, 0.52, 0.4, 0.4], 'String', ['vx=',num2str(vx), ...
            char(10) 'vy = ',num2str(vy), ...
		    char(10) 'T0 =',num2str(T0),' [K]' ...lin
		    char(10) 'q = ',num2str(q),' [W]' ...
		    char(10) 'k = ',num2str(k),' [W/m/K]' ...
		    char(10) 'cp = ',num2str(cp),' [J/kg/K]' ...
		    char(10) 'rho = ',num2str(rho),' [kg/m^3]'],'LineStyle','none','color','g');

for n = 1:1:Data_Points
    T_Array_Trial = Trial_Temps(:,((n-1)*nx)+1:n*nx);     %Selects the Data for given time
    % Creates a contour plot for the data
    [Cont,h]=contourf(x,y,T_Array_Trial,[T0:10:Max_Temp*(Temp_sensitivity/100)],'LineStyle','none');
    cmap = jet(256); % Sets the color profile of the contour plot
    colormap(cmap);
    view(2);
    colorbar
    axis tight 
    caxis([T0,Max_Temp*(Temp_sensitivity/100)]);

    xlabel('x [m]') % Updates labels 
    ylabel('y [m]')
    title(['Rosenthals solution T: for t = ',num2str(Time(n))])
    xlim([-((x_length/2)+Window_Pos(1))*Scale, ((x_length/2)+Window_Pos(1))*Scale])
    ylim([-((y_length/2)+Window_Pos(2))*Scale, ((y_length/2)+Window_Pos(2))*Scale])
    pause(0.0001) % Ensures that the plots show up
end       