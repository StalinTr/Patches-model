
function [N, P, T] = patches(domain, tspan, params, eddy, iCondN, iCondP, iCondType)
% PATCHES   Simulates pytoplankton patches in a 2D plane
% 
%   IN:
%       domain:     Space domain. A 6x1 vector containing:
%                     - x-coordinate of the lower left corner (m)
%                     - y-coordinate of the lower left corner (m)
%                     - x-coordinate of the upper right corner (m)
%                     - y-coordinate of the upper right corner (m)
%                     - Number of cells in the x direction
%                     - Number of cells in the y direction
% 
%       tspan:      Time span. A 1x3 vector containing:
%                     - Start time (days)
%                     - Step size (days)
%                     - End time (days)
% 
%       params:     Parameter struct containing the physical parameters as
%                   given in the following table:
%
%                   VARIABLE  | UNITS         | DESCRIPTION
%                  -----------+---------------+-------------------------------------------------------
%                    IN       |  g·m-2·day-1  |  Input rate of nutrient
%                    mu       |  day-1        |  Maximum growth rate of phytoplankton
%                    k        |               |  Nutrient content in phytoplankton
%                    mN       |  day-1        |  Loss rate of nutrient
%                    fP       |  g·m-2·day-1  |  Maximum feeding rate of zooplankton on phytoplankton
%                    HN       |  g·m-2        |  Half-saturation constant of nutrient
%                    HP       |  g·m-2        |  Half-saturation constant of phytoplankton
%                    d        |  ???          |  Diffusion coefficient
%                    vmax     |  m/s          |  Maximum fluid velocity
%                    A        |  %            |  Amplitude as a percentage of distance from
%                             |               |  stationary point to origin (if iCondType == 'scaled')
%   
%       eddy:       Eddy currents. Struct containing the following fields as Nx1
%                   vectors, where N is the number of eddies of the simulation.
%                   NOTE: Mirror eddies will be created internally to enforce 
%                   periodic boundary conditions.
%                   
%                   VARIABLE  | UNITS         | DESCRIPTION
%                  -----------+---------------+--------------------------------------------------------
%                    x        |  m            |  x-coordinate of the center of the eddies
%                    y        |  m            |  y-coordinate of the center of the eddies
%                    r        |  m            |  Radius of the eddies
%                    s        |               |  Rotation directon (1: clockwise, -1: counter-clockwise)
%
%       iCondN:     Initial conditions. Initial concentration of nutrients (in g·m3)
%                   as a NxM matrix where N is the number of cells in the y direction
%                   and M the number of cells in the x direction. Initial conditions
%                   must be centered around a instable stationary point in the
%                   dimensionless space, and must have periodic boundaries.
% 
%       iCondP:     Initial conditions. Initial concentration of phytoplankton (in g·m3)
%                   as a NxM matrix where N is the number of cells in the y direction
%                   and M the number of cells in the x direction. Initial conditions
%                   must be centered around a instable stationary point in the
%                   dimensionless space, and must have periodic boundaries.
%
%       iCondType:  Initial conditions type. This is an optional parameter. If it's value
%                   is 'adjust' then the initial conditions will be linearly mapped so
%                   that their mean is located at a valid stationary point and there are
%                   no negative values. Note the initial conditions are only scaled if
%                   negative values are found after translation. This option aims to
%                   help in choosing valid initial conditions.
% 
%   OUT:
%       N:          Nutrient concentration in a NxMxL matrix, where N is the 
%                   number of cells in the y direction, M the number of cells in
%                   the x direction and L the number of time steps.
% 
%       P:          Pytoplankton concentration in a NxMxL matrix, where N is the 
%                   number of cells in the y direction, M the number of cells in
%                   the x direction and L the number of time steps.
% 


    %==========================================================================
    % PARAMETERS
    %==========================================================================

    dt = tspan(2);      % Time step (0.12-0.36)

    % Create time vector
    T = tspan(1):dt:tspan(3);   % Actual time
    t = params.mu * T;          % Dimensionless time
    
    in = params.IN / (params.mu * params.HN);   % Dimensionless input rate of nutrients 0.9, 0.32, 0.38
    a  = params.k * params.HP / (params.HN);    % Dimensionsless nutrient content in phytoplankton 0.8
    mn = params.mN / params.mu;                 % Dimensionless loss rate of nutrient 0.03
    fp = params.fP / (params.mu * params.HP);   % Dimensionless max feeding rate of zoo on phyto 0.9

    %==========================================================================
    % SETUP
    %==========================================================================

    % Size of the domain
    Lx = domain(3) - domain(1);
    Ly = domain(4) - domain(2);
    
    % Space steps
    dx = Lx / domain(5);
    dy = Ly / domain(6);
    
    % Create grid
    [x,y] = meshgrid( linspace( domain(1), domain(3), domain(5) ),...
                      linspace( domain(2), domain(4), domain(6) ) );
    
    % Add mirror eddies so that the stream field has periodic boundaries
    offx = Lx + dx;         % Offset of mirrors in the x direction
    offy = Ly + dy;         % Offset of mirrors in the y direction
    eddy.x = [eddy.x; offx+eddy.x; -offx+eddy.x; eddy.x; eddy.x; offx+eddy.x; offx+eddy.x; -offx+eddy.x; -offx+eddy.x;];
    eddy.y = [eddy.y; eddy.y; eddy.y; offy+eddy.y; -offy+eddy.y; offy+eddy.y; -offy+eddy.y; offy+eddy.y; -offy+eddy.y];
    eddy.s = repmat(eddy.s,[9,1]);
    eddy.r = repmat(eddy.r,[9,1]);
    
    % Get preliminary stream field and velocities
    psi = psiFun(x,y);
    vx = -Dy(psi);
    vy = Dx(psi);
    v = sqrt(vx.^2+vy.^2);

    % Adjust s parameter so that the max velocity is vmax
    mult = params.vmax / max(v,[],'all');
    vx = vx*mult;
    vy = vy*mult;
    
    % Calculate velocity derivatives to use in calculations
    dv = Dx(vx) + Dy(vy);

    % Initialize variables
    n = zeros([size(x),length(t)]); % Dimensionless nutrient concentration
    p = zeros([size(x),length(t)]); % Dimensionless phytoplankton concentration    

    % Initial conditions
    % iCondType is given and is 'adjust'
    if nargin > 6 && strcmpi(iCondType, 'adjust')
        % Adjust initial conditions and assign to dimensionless conditions
        [iCondn, iCondp] = adjustIC(iCondN, iCondP);
        n(:,:,1) = iCondn;
        p(:,:,1) = iCondp;
        
    % iCondType is not given or is not 'adjust'
    else
        % Assign initial dimensionless conditions
        n(:,:,1) = iCondN / params.HN;
        p(:,:,1) = iCondP / params.HP;
    end
    
    % Error detection. Concentrations should not go below zero. Check the given parameters.
    if any(n(:,:,1) < 0, 'all')
        warning('The initial conditions for N (iCondN) has negative concentrations. Consider adjusting the input parameters.');
    end
    if any(p(:,:,1) < 0, 'all')
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
        [nTemp, pTemp] = RK4(n(:,:,tau-1), p(:,:,tau-1));
        n(:,:,tau) = nTemp;
        p(:,:,tau) = pTemp;

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
    % Computes the first-order central finite difference on x
        dA = ([A(:,2:end),A(:,1)]-[A(:,end),A(:,1:end-1)])/(2*dx);
    end

    function [ddA] = DDx(A)
    % Computes the second-order central finite difference on x
        ddA = ([A(:,2:end),A(:,1)]-2*A+[A(:,end),A(:,1:end-1)])/dx^2;
    end

    function [dA] = Dy(A)
    % Computes the first-order central finite difference on y
        dA = ([A(2:end,:);A(1,:)]-[A(end,:);A(1:end-1,:)])/(2*dy);
    end

    function [ddA] = DDy(A)
    % Computes the second-order central finite difference on y
        ddA = ([A(2:end,:);A(1,:)]-2*A+[A(end,:);A(1:end-1,:)])/dy^2;
    end

    function [dn] = DnFun(n,p)
    %     dn = d*(DDx(n)+DDy(n)) - n.*dv - (vx.*Dx(n)+vy.*Dy(n)) + in - a*n./(1+n).*p - mn*n;
        dn = 0;
        dn = dn - n.*dv - (vx.*Dx(n)+vy.*Dy(n));                       % Advection
        dn = dn + params.d*(DDx(n)+DDy(n));                            % Diffusion
        dn = dn + in - a*n./(1+n).*p - mn*n;                           % Reaction
    end

    function [dp] = DpFun(n,p)
    %     dp = d*(DDx(p)+DDy(p)) - p.*dv - (vx.*Dx(p)+vy.*Dy(p)) + n./(1+n).*p - fp*p./(1+p);
        dp = 0;
        dp = dp - p.*dv - (vx.*Dx(p)+vy.*Dy(p));                       % Advection
        dp = dp + params.d*(DDx(p)+DDy(p));                            % Diffusion
        dp = dp + n./(1+n).*p - fp*p./(1+p);                           % Reaction
    end

    function [psi] = psiFun(x,y)
    % Computes the psi function at the positions given by matrices x and y

        % Check that sizes match
        if any(size(x) ~= size(y))
            error('x and y must have the same size');
        end

        s = 1; % Scaling factor, adjusted in script
        psi = zeros(size(x));
        for i = 1:size(x,1)
            for j = 1:size(x,2)
                psi(i,j) = s*sum(eddy.s.*exp(-((x(i,j)-eddy.x).^2+(y(i,j)-eddy.y).^2)./eddy.r.^2));
            end
        end
    end

    function [n,p] = RK4(n,p)
    % Runge-Kutta 4 for f(x,y)

        k1 = dt*DnFun(n,p);
        m1 = dt*DpFun(n,p);

        k2 = dt*DnFun(n+k1/2, p+m1/2);
        m2 = dt*DpFun(n+k1/2, p+m1/2);

        k3 = dt*DnFun(n+k2/2, p+m2/2);
        m3 = dt*DpFun(n+k2/2, p+m2/2);  

        k4 = dt*DnFun(n+k3, p+m3);
        m4 = dt*DpFun(n+k3, p+m3);

        n = n + (k1+(k2*2)+(k3*2)+k4)/6;
        p = p + (m1+(m2*2)+(m3*2)+m4)/6;

    end

end