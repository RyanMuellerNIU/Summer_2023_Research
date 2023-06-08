% Rosenthal's solution: Written By Jakub Mikula
% Input parameters:
	T0 = 300; 			% ambient temperature [K]
	q = 840; 			% laser power [W]
	k = 35; 			% conductivity coefficient [W/m/K]
	v = 1000/1000/60; 	% laser speed [m/s]
	cp = 800; 			% specific heat [J/kg/K]
	rho = 7600.0; 		% density [kg/m^3]
	nx = 1000;   % number of grid points
	ny = 1000;   % number of grid points 

	[x,y] = meshgrid(linspace(-1,1,nx)/5,linspace(-1,1,nx)/5); % [m]

			for i=1:nx
			for j=1:ny
				a = k/rho/cp;
				w = x(i,j);
				R = sqrt(w^2+y(i,j)^2);
				T(i,j) = T0+q/2/pi/k/R*exp(-v*(w+R)/2/a);
			end
            end
            T_Max = max(max(T));
            figure(2)
			[C,h]=contourf(x,y,T,[T0:10:1500],'LineStyle','none');
			view(2);
			caxis([T0,1500]);
			clabel(C,h,[T0:10:400],'FontSize',6,'Color','white');

	annotation('textbox', [0.13, 0.52, 0.4, 0.4], 'String', ['v=',num2str(v), ...
		char(10) 'T0 =',num2str(T0),' [K]' ...lin
		char(10) 'q = ',num2str(q),' [W]' ...
		char(10) 'k = ',num2str(k),' [W/m/K]' ...
		char(10) 'cp = ',num2str(cp),' [J/kg/K]' ...
		char(10) 'rho = ',num2str(rho),' [kg/m^3]'],'LineStyle','none','color','w');

	colorbar
	axis tight 

	xlabel('x [m]')
	ylabel('y [m]')
	title('Rosenthals solution T')

	xlim([-0.05, 0.05])
	ylim([-0.05, 0.05])