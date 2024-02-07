clear; clc

Imax        = 800;      % | Wm-2      | Maximum irradiance
IOffset     = 200;      % | Wm-2      | Irradiance generation offset
Tday        = 24;       % | h         | Day duration
tNoon       = 12;       % | h         | Time to noon (highest irradiance)
Ki          = 0.5;      % | m-1       | Light atenuation coefficient 0.8


% Light irradiance at surface level
    Is    = @(t) max( 0, (Imax-IOffset)*sin(2*pi*(t-tNoon)/Tday+pi/2)+IOffset );

% Light transmission within the water
getIhandle = @(Is, z, Ki) getI(Is, z, Ki);

z = -10:0.1:0;
t = 0:0.1:24;

figure(1); clf
plotPreLightIrradiance(t, Ki, Imax, Is, getIhandle)

function [I] = getI(Is, z, Ki)
    % Computes the light irradiance that reaches each cell given a surface
    % irradiance and the position of the particles. It uses the
    % Beer-Lambert attenuation model.
    %
    % INPUT     Is      (1x1) Light irradiance at the surface.
    %           z       (Nx1) Vertical position of the cells.
    %
    % OUTPUT    I       (Nx1) Light irradiance at the position of each cell

    I = Is*exp(Ki*z);

end

function [] = plotPreLightIrradiance(tVec, Ki, Imax, Is, getI)
    % Plots 3 subplots:
    %   1.  The light irradiance that reaches the surface of the water over
    %       time.
    %   2.  The light irradiance that reaches different depths over time.
    %   3.  The maximum light irradiance that is reached at different
    %       depths.
    %
    % INPUT     tVec    (1xN) Time vector (h).
    %           Ki      (1x1) Light atenuation coefficient (m-1).
    %           Ilim    (1x1) Irradiance thresshold for light/dark conditions
    %           Imax    (1x1) Maximum irradiance (Wm-2)
    %           Is      (func. handle) Light intensity that reaches the
    %                   surface (Nx1). It is a function of time (Nx1).
    %           getI    (func. handle) Light intensity that reaches
    %                   different depths. See getI in simu.m
    
    % Plot parameters
    propX  = 0.8;   propY  = 0.3;
    offX   = 0.025; offY   = 0.1;
    margX  = 0.05;  margY  = 0.15;
    interX = 0.02;  interY = 0.02;
    fontSize = 10;
    tLim = [tVec(1), tVec(end)];

    % Compute light irradiance at depths and times
    tN = length(tVec);
    dt = tVec(2)-tVec(1);
    zVecPlot = -15:0.01:0;
    II = zeros(tN,length(zVecPlot));
    for n = 1:tN
        t = tVec(n);
        II(n,:) = getI(Is(t), zVecPlot, Ki);
    end 

    % Plot 1: Surface light irradiance over time
    h1 = subplot(3,1,1);
        plot(tVec, Is(tVec));
        xlim(tLim);
        ylim([-99, 890]);
        ylabel('I_s (W/m^2)','FontSize',fontSize);
        set(gca,'xticklabel',[]);

    % Plot 2: Light irrandiance at depths and over time
    h2 = subplot(3,1,2);
        sh = pcolor(tVec, zVecPlot, II'); hold on   % Color light intensity
        xlim(tLim);
        c = colorbar('southoutside');
        c.Label.String = 'I (W/m^2)';
        set(sh,'EdgeColor', 'none' );
        xlabel('t (h)','FontSize',fontSize);
        ylabel('z (m)','FontSize',fontSize);
   
    % Plot 3: Maximum light irradiance that reached each depth
    h3 = subplot(3,1,3);
        II= getI(Imax, zVecPlot, Ki);
        zVecPlot = [zVecPlot, dt:dt:4.5];%dt:dt:3.05];
        II = [II, Imax*ones(1, length(zVecPlot)-length(II))];
        plot(II, zVecPlot); 
        xlim([-99,890]); ylim([min(zVecPlot), max(zVecPlot)]); yticks(-10:2:2); yticklabels([]);
        xlabel('max(I) (W/m^2)','FontSize',fontSize);%ylabel('z (m)','FontSize',fontSize);

    % Reposition axis

        % Plot 1
        hpos = [offX+margX, offY+1-propY+1/2*interY, propX-margX-1/2*interX, propY-margY-1/2*interY]; 
        set(h1, 'position', hpos);

        % Plot 2
        hpos = [offX+margX, offY+margY, propX-margX-1/2*interX, 1-propY-margY-1/2*interY];
        set(h2, 'position', hpos);
        
        % Plot3
        hpos = [offX+propX+1/2*interX, offY+margY, 1-propX-margX-1/2*interX, 1-2*margY];
        set(h3, 'position', hpos);

end