
function [N, P, T] = patches3D(domain, tspan, params, eddy, iCondN, iCondP, iCondType)
	% PATCHES3D   Simulates pytoplankton patches in a 3D space
	
		%==========================================================================
		% PARAMETERS
		%==========================================================================
	
		dt = tspan(2);      % Time step (0.12-0.36)
	
		% Create time vector
		T = tspan(1):dt:tspan(3);   % Actual time
		t = params.mu * T;          % Dimensionless time
		
		in = params.IN / (params.mu * params.HN);   % Dimensionless input rate of nutrients
		a  = params.k * params.HP / (params.HN);    % Dimensionsless nutrient content in phytoplankton
		mn = params.mN / params.mu;                 % Dimensionless loss rate of nutrient
		fp = params.fP / (params.mu * params.HP);   % Dimensionless max feeding rate of zoo on phyto
	
		%==========================================================================
		% SETUP
		%==========================================================================
	
		% Size of the domain
		Lx = domain(4) - domain(1);
		Ly = domain(5) - domain(2);
		Lz = domain(6) - domain(3);
		
		% Space steps
		dx = Lx / domain(7);
		dy = Ly / domain(8);
		dz = Lz / domain(9);
		
		% Create grid
		[x,y,z] = meshgrid( linspace( domain(1), domain(4), domain(7) ),...
						  linspace( domain(2), domain(5), domain(8) ), ...
						  linspace( domain(3), domain(6), domain(9) ) );
		
		% Add mirror eddies so that the stream field has periodic boundaries
		offx = Lx + dx;         % Offset of mirrors in the x direction
		offy = Ly + dy;         % Offset of mirrors in the y direction
		offz = Lz + dz;         % Offset of mirrors in the z direction
		eddy.x = [eddy.x; offx+eddy.x; -offx+eddy.x; eddy.x; eddy.x; offx+eddy.x; offx+eddy.x; -offx+eddy.x; -offx+eddy.x;];
		eddy.y = [eddy.y; eddy.y; eddy.y; offy+eddy.y; -offy+eddy.y; offy+eddy.y; -offy+eddy.y; offy+eddy.y; -offy+eddy.y];
		eddy.z = [eddy.z; eddy.z; eddy.z; eddy.z; eddy.z; eddy.z; eddy.z; eddy.z; eddy.z];
		eddy.s = repmat(eddy.s,[9,1]);
		eddy.r = repmat(eddy.r,[9,1]);
		
		% Get preliminary stream field and velocities
		psi = psiFun(x,y,z);
		vx = -Dy(psi);
		vy = Dx(psi);
		vz = zeros(size(psi));
		v = sqrt(vx.^2 + vy.^2 );
	
		% Adjust s parameter so that the max velocity is vmax
		mult = params.vmax / max(v,[],'all');
		vx = vx * mult;
		vy = vy * mult;
		vz = vz * mult;
		
		% Calculate velocity derivatives to use in calculations
		dv = Dx(vx) + Dy(vy); %+ Dz(vz);
	
		% Initialize variables
		n = zeros([size(x),length(t)]); % Dimensionless nutrient concentration
		p = zeros([size(x),length(t)]); % Dimensionless phytoplankton concentration    
	
		% Initial conditions
		% iCondType is given and is 'adjust'
		if nargin > 6 && strcmpi(iCondType, 'adjust')
			% Adjust initial conditions and assign to dimensionless conditions
			[iCondn, iCondp] = adjustIC(iCondN, iCondP);
			n(:,:,:,1) = iCondn;
			p(:,:,:,1) = iCondp;
			
		% iCondType is not given or is not 'adjust'
		else
			% Assign initial dimensionless conditions
			n(:,:,:,1) = iCondN / params.HN;
			p(:,:,:,1) = iCondP / params.HP;
		end
		
		% Error detection. Concentrations should not go below zero. Check the given parameters.
		if any(n(:,:,:,1) < 0, 'all')
			warning('The initial conditions for N (iCondN) has negative concentrations. Consider adjusting the input parameters.');
		end
		if any(p(:,:,:,1) < 0, 'all')
			warning('The initial conditions for P (iCondP) has negative concentrations. Consider adjusting the input parameters.');
		end
		
		%==========================================================================
		% CALCULATIONS
		%==========================================================================
		
		fprintf('Time steps:   0%%');
		tN = length(t);
			
		for tau = 2:tN
			
			if mod(tau, round(tN/100))==0
				fprintf('\b\b\b\b%3i%%', round(tau/tN*100));
			end
				
			% Runge-Kutta 4
			[nTemp, pTemp] = RK4(n(:,:,:,tau-1), p(:,:,:,tau-1));
			n(:,:,:,tau) = nTemp;
			p(:,:,:,tau) = pTemp;
	
		end
		
		fprintf('\b\b\b\b%3i%%', 100);
	
		% Error detection. Concentrations should not go below zero. Check the given parameters.
		if any(n < 0, 'all')
			warning('The solution for N has negative concentrations. Consider adjusting the input parameters.');
		end
		if any(p < 0, 'all')
			warning('The solution for P has negative concentrations. Consider adjusting the input parameters.');
		end

        % Transform dimensionless representation to physical
        N = n * params.HN;
        P = p * params.HP;

    
        fprintf('\nDone!\n')
	
		%==========================================================================
		% FUNCTIONS
		%==========================================================================

    function [iCondn, iCondp] = adjustIC(iCondN, iCondP)
    % Adjusts initial conditions so that the concentration distribution in
    % the dimensionless space is centered on a valid stationary point and
    % does not have any negative values.
    
        % Compute stationary point in dimensionless space (n1, p1)
        %   There are 3 SP: (n0, 0), (n1, p1) and (n2, p2).
        %   The first is no good since there would be no phytoplankton.
        %   The other 2 are solutions to
        %       n^2+c2n+c3=0,
        %       p=fp/n+fp-1,
        %   with c1, c2,c3 given below. Also, n and p must be greater than 0 to
        %   seed the initial condition.
        c1 = 1; c2 = (a*fp-in-a+mn)/mn; c3 = (a*fp-in)/mn;
        n1 = (-c2+sqrt(c2^2-4*c1*c3))/(2*c1);
        p1 = fp/n1+fp-1;
        % The first stationary point was negative in at least one of the axes
        if n1 < 0 || p1 < 0
            err(1,:) = [n1,p1];
            n1 = (-c2-sqrt(c2^2-4*c1*c3))/(2*c1);
            p1 = fp/n1+fp-1;
        end
        % Both stationary points were negative in at least one of the axes
        if n1 < 0 || p1 < 0
            err(2,:) = [n1,p1];
            disp(err)
            error('The given parameters do not have a valid stationary point')
        end
        
        % Transform to dimensionless space
        iCondn = iCondN / params.HN;
        iCondp = iCondP / params.HP;
        
        % Translate so that the mean is centered at the stationary point
        iCondn = iCondn + n1 - mean(iCondn, 'all');
        iCondp = iCondp + p1 - mean(iCondp, 'all');
        
        % Scale range keeping the mean if there are negative values
        if any(iCondn <= 0, 'all')
            iCondn = (iCondn - n1)/(n1-min(iCondn, [], 'all'))*(params.A*n1)+n1;    
        end
        if any(iCondp <= 0, 'all')
            iCondp = (iCondp - p1)/(p1-min(iCondp, [], 'all'))*(params.A*p1)+p1;    
        end
        
        % Display adjustment in physical space
        iCondN = iCondn * params.HN; iCondP = iCondp * params.HP;
        fprintf('iCondN adjusted to [%f, %f] with mean %f\n', min(iCondN,[],'all'), max(iCondN,[],'all'), mean(iCondN,'all'));
        fprintf('iCondP adjusted to [%f, %f] with mean %f\n', min(iCondP,[],'all'), max(iCondP,[],'all'), mean(iCondP,'all'));
        
        % Error checking
        if abs(mean(iCondn,'all') - n1) > 1e-10
            warning('Something went wrong in adjustment. Mean of iCondn does not match stationary point n1 (%f vs. f%)', mean(iCondn,'all'), n1)
        end
        if abs(mean(iCondp,'all') - p1) > 1e-10
            warning('Something went wrong in adjustment. Mean of iCondp does not match stationary point p1 (%f vs. f%)', mean(iCondp,'all'), p1)
        end
        
    end
	
		function [dA] = Dx(A)
			% Computes the derivative of A in the x direction
			dA = (A([2:end,1],:,:) - A([end,1:end-1],:,:)) / (2*dx);
		end

		function [dA] = Dy(A)
			% Computes the derivative of A in the y direction
			dA = (A(:,[2:end,1],:) - A(:,[end,1:end-1],:)) / (2*dy);
		end

		function [dA] = Dz(A)
			% Computes the derivative of A in the z direction
			dA = (A(:,:,[2:end,1]) - A(:,:,[end,1:end-1])) / (2*dz);
		end

		function [dA] = DDx(A)
			% Computes the second derivative of A in the x direction
			dA = (A([2:end,1],:,:) - 2*A + A([end,1:end-1],:,:)) / dx^2;
		end

		function [dA] = DDy(A)
			% Computes the second derivative of A in the y direction
			dA = (A(:,[2:end,1],:) - 2*A + A(:,[end,1:end-1],:)) / dy^2;
		end

		function [dA] = DDz(A)
			% Computes the second derivative of A in the z direction
			dA = (A(:,:,[2:end,1]) - 2*A + A(:,:,[end,1:end-1])) / dz^2;
		end

		function [dn] = DnFun(n,p)
			% dn = d*(DDx(n)+DDy(n)+DDz(n)) - n.*dv - (vx.*Dx(n)+vy.*Dy(n) +vz.*Dz(n)) + in - a*n./(1+n).*p - mn*n;
			dn = 0;
			dn = dn - n.*dv - (vx.*Dx(n)+vy.*Dy(n) +vz.*Dz(n));	% Advection
			dn = dn + params.d*(DDx(n)+DDy(n)+DDz(n));		% Diffusion
			dn = dn + in - mn*n - a*n./(1+n).*p; % Reaction
		end

		function [dp] = DpFun(n,p)
			% dp = d*(DDx(p)+DDy(p)+DDz(p)) - p.*dv -(vx.*Dx(p)+vy.*Dy(p) +vz.*Dz(p)) + n./(1+n).*p - fp*p./(1+p);
			dp = 0;
			dp = dp - p.*dv - (vx.*Dx(p)+vy.*Dy(p) +vz.*Dz(p));	% Advection
			dp = dp + params.d*(DDx(p)+DDy(p)+DDz(p));	% Diffusion
			dp = dp + n./(1+n).*p - fp*p./(1+p);	% Reaction
		end
	
		function [psi] = psiFun(x,y,z)
		% Computes the psi function at the positions given by matrices x,y and z
			
			% Check that sizes match
			if any(size(x) ~= size(y)) || any( size(x) ~= size(z))
				error('x, y or z must have the same size');
			end
			
			s = 1; % Scaling factor, adjusted in script
			psi = zeros(size(x));
			for i = 1:size(x,1)
				for j = 1:size(x,2)
					for k = 1:size(x,3)
						psi(i,j,k) = s*sum(eddy.s.*exp(-((x(i,j,k)-eddy.x).^2+(y(i,j,k)-eddy.y).^2 + (z(i,j,k) -eddy.z).^2)./eddy.r.^2));
					end	
				end
			end
		end

		function [n, p] = RK4(n, p)
			% Runge-Kutta 4
			k1n = dt*DnFun(n,p);
			k1p = dt*DpFun(n,p);

			k2n = DnFun(n+k1n/2, p+ k1p/2);
			k2p = DpFun(n+k1n/2, p+ k1p/2);

			k3n = DnFun(n+k2n/2, p+k2p/2);
			k3p = DpFun(n+k2n/2, p+k2p/2);

			k4n = DnFun(n+k3n, p+k3p);
			k4p = DpFun(n+k3n, p+k3p);

			n = n + (k1n+2*k2n+2*k3n+k4n)/6;
			p = p + (k1p+2*k2p+2*k3p+k4p)/6;
		end

	
	end
	