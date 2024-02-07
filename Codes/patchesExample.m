clear; clc

vidName = 'video.mp4';

domain = [-400,-200,400,200,200,100];
Lx = domain(3)-domain(1);
Ly = domain(4)-domain(2);

tSpan = [0, 0.2, 185];

params.IN   = 0.00225;% 0
params.mu   = 0.5;
params.k    = 0.01;
params.mN   = 0.015;
params.fP   = 1.8; % 2.0
params.HN   = 0.005;
params.HP   = 4.0;
params.d    = 0.04;
params.vmax = 0.3;
params.A   = 0.5;

nEddy = 100;
eddy.x = Lx*rand(nEddy, 1)+domain(1);       % x_i
eddy.y = Ly*rand(nEddy, 1)+domain(2);       % y_i
eddy.s = 2*(round(rand(nEddy, 1))-0.5);     % sigm_i
eddy.r = abs(normrnd(20,5,[nEddy,1]));      % r_i

[x,y] = meshgrid( linspace( domain(1), domain(3), domain(5) ),...
                      linspace( domain(2), domain(4), domain(6) ) );

% Initial conditions
% Must be near (n1, p1) and should have periodic boundary conditions
iCondN = (sin(pi/(Lx/2)*(x-100))).*(sin(pi/(Ly/2)*(y-0)));
iCondP = sqrt(x.^2+y.^2);

%iCondN = zeros(size(iCondN));
%iCondP = iCondP - mean(iCondP,'all');
%iCondP = iCondP / max(iCondP,[],'all')*0.4563/2;
%iCondP = iCondP + 0.484514;

iCondN = 1+params.A*sin(pi/Lx*2*x);
iCondP = 1+params.A*sin(pi/Ly*2*y);

%xc = 0; yc = 0; sigma = 5000;
%iCondP = (exp(-((x-xc).^2 + (y-yc).^2)./(2*sigma)));
%fprintf('iCondP is in [%f, %f] with mean %f\n', min(iCondP,[],'all'), max(iCondP,[],'all'), mean(iCondP,'all'));

[N, P, T] = patches(domain, tSpan, params, eddy, iCondN, iCondP,'adjust'); % 'adjust'

%==========================================================================
% Plots
%==========================================================================

% Color map
cmapC = [0.023529411764706   0.070588235294118   0.388235294117647;
         0.509803921568627   1.000000000000000   0.509803921568627;
         0.933333333333333   1.000000000000000   0.349019607843137];
cmapP = [0,0.9,1];
cmap = interp1(cmapP,cmapC,(1:255)/255);

% Eddy currents
figure(1); clf
    scatter([], [], 'r'); hold on   % Used in legend
    scatter([], [], 'b');
    viscircles([eddy.x(eddy.s == 1), eddy.y(eddy.s == 1)], eddy.r(eddy.s == 1), 'Color', 'r');
    viscircles([eddy.x(eddy.s == -1), eddy.y(eddy.s == -1)], eddy.r(eddy.s == -1), 'Color', 'b');
    rectangle('Position',[domain([1,2]), domain([3,4])-domain([1,2])]); hold off
    title('Eddy currents');
    legend('Clockwise', 'Counter-clockwise');
    axis equal tight

% Initial and final state
figure(2); clf
    subplot(2,2,1)
        h = pcolor(x, y, N(:,:,1));     set(h, 'EdgeColor', 'none'); colorbar; colormap(cmap);
        title(['N (t = ', num2str(tSpan(1),'%.2f'), ' days)']); axis equal tight
    subplot(2,2,2)
        h = pcolor(x, y, N(:,:,end));   set(h, 'EdgeColor', 'none'); colorbar; colormap(cmap);
        title(['N (t = ', num2str(tSpan(3),'%.2f'), ' days)']); axis equal tight
    subplot(2,2,3)
        h = pcolor(x, y, P(:,:,1));     set(h, 'EdgeColor', 'none'); colorbar; colormap(cmap);
        title(['P (t = ', num2str(tSpan(1),'%.2f'), ' days)']); axis equal tight
    subplot(2,2,4)
        h = pcolor(x, y, P(:,:,end));   set(h, 'EdgeColor', 'none'); colorbar; colormap(cmap);
        title(['P (t = ', num2str(tSpan(3),'%.2f'), ' days)']); axis equal tight

% Periodicity check
figure(3); clf
    subplot(2,2,1); 
        imagesc(repmat(N(:,:,1),[3,3]));    title(['N (t = ', num2str(tSpan(1),'%.2f'), ' days)']);
        set(gca,'YDir','normal'); colormap(cmap); axis equal tight off
    subplot(2,2,2);
        imagesc(repmat(N(:,:,end),[3,3]));  title(['N (t = ', num2str(tSpan(3),'%.2f'), ' days)']);
        set(gca,'YDir','normal'); colormap(cmap); axis equal tight off
    subplot(2,2,3);
        imagesc(repmat(P(:,:,1),[3,3]));    title(['P (t = ', num2str(tSpan(1),'%.2f'), ' days)']);
        set(gca,'YDir','normal'); colormap(cmap); axis equal tight off
    subplot(2,2,4);
        imagesc(repmat(P(:,:,end),[3,3]));  title(['P (t = ', num2str(tSpan(3),'%.2f'), ' days)']);
        set(gca,'YDir','normal'); colormap(cmap); axis equal tight off
    sgtitle('Periodicity check')

% Record video
    vidFrames = round(linspace(1, length(T), 200));
    recordVideo(vidName, vidFrames, N, P, T, x, y, cmap)


function [] = recordVideo(vidName, vidFrames, N, P, T, x, y, cmap)
% Records a video of the results skipping vidSkip frames each time

    % Decide wether to overwrite output file if already exists
    if exist(vidName, 'file')==2
        ufig = uifigure;
        msg = ['Runing the script will overwrite a file: ',vidName];
        tit = 'Confirm Save';
        selection = uiconfirm(ufig,msg,tit, ...
               'Options',{'Overwrite','Save as new','Cancel'}, ...
               'DefaultOption',2);
        if strcmp(selection, 'Overwrite')
            delete(vidName);
        elseif strcmp(selection, 'Save as new')
            i = 1;
            vidName = [vidName(1:end-4), '(', num2str(i), ').mp4'];
            while exist(vidName, 'file')==2
                nameMatches = regexp(vidName,'(\d)+');
                i = i+1;
                vidName = [vidName(1:nameMatches-2), '(', num2str(i), ').mp4'];
            end
        elseif strcmp(selection, 'Cancel')
            close(ufig)
            return
        end
        close(ufig) 
    end
    
    % Open video file
    vid = VideoWriter(vidName,'MPEG-4');
    vid.Quality = 100;
    open(vid);
      
    % Write video
    figure(99); clf; hold off
    fprintf('Recording video:   0%%');
    for i = vidFrames

        if mod(i, round(length(vidFrames)/100))==0
            fprintf('\b\b\b\b%3i%%', round(i/length(vidFrames)*100));
        end
        
        subplot(2,1,1);
            h = pcolor(x, y, N(:,:,i));
            set(h, 'EdgeColor', 'none'); 
            title('N'); 
            colorbar;
            colormap(cmap);
            axis equal tight
        subplot(2,1,2);
            h = pcolor(x, y, P(:,:,i));
            set(h, 'EdgeColor', 'none');
            title('P');
            colorbar;
            colormap(cmap);
            axis equal tight
        
        sgtitle(['t = ', num2str(T(i), '%.2f'), ' days']);
        %set(gca,'xtick',[])
        %set(gca,'ytick',[])
        drawnow
        frame = getframe(gcf);
        writeVideo(vid,frame);
    end
    
    fprintf('\b\b\b\b%3i%%', 100);

    % Close video file
    close(vid);

    fprintf('\nDone!\n')
end
