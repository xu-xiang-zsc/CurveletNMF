function [U,P] = OSP(HIM,p)

% OSP algorithm for endmember extraction.
% ------------------------------------------------------------------------------
% Input:   HIM : hyperspectral image cube [nrows x ncols x nchannels]
%          p   : desired number of endmembers to be extracted
%
% Output:  U   : set of extracted endmembers [nchannels x p] 
%          P   : spatial coordinates of the extracted endmembers (positions
%                rows x cols)
% 
% Copyright (2007) GRNPS group @ University of Extremadura, Spain. 


disp(' ====== Start OSP run ======');

% get current CPU time
t1 = cputime;

% get image size
[ns,nl,nb] = size(HIM);

if nb == 1
    Uin_row(1,:,:) = HIM';
    HIM = Uin_row;
    [ns,nl,nb] = size(HIM);
end

% visualize image
imagesc(mean(HIM,3)); colormap(gray); 
set(gca,'DefaultTextColor','black','xtick',[],'ytick',[],'dataaspectratio',[1 1 1]);
% po1 = get(gca,'position');

% Calculate the pixel (vector) with major intensity in the image
max = 0;
for i = 1 : ns
    for j = 1 : nl
        r = squeeze(HIM(i,j,:));
        bright = r'*r;
        if bright > max
            max = bright;
            posx = i;
            posy = j;
        end
    end
end

% The pixel with more intensity is the initial pixel of the process
t0 = squeeze(HIM(posx,posy,:));

% Generate the identity matrix.
I = eye(nb,nb);

% Initialization of the pure pixels matrix
U = [];
U = [U t0];

% Initialization of the positions matrix
P = zeros(p, 2);
P(1,1)=posx; P(1,2)=posy;

fprintf(1, '. found pixel @ coordinates x=%5d & y=%5d \n',posx,posy);

% Visualization of the position of the first chosen pixel
drawnow;
text(posy,posx,'o','Margin',1,'HorizontalAlignment','center','FontSize',22,'FontWeight','light','FontName','Garamond','Color','yellow');
% text(posy,posx,'o','Color','yellow');


%% Algoritmo OSP
for i = 1:p-1
    UC = U(:,1:i);

    % Calculate the orthogonal projection with respect to the pixels at present chosen. 
    % This part can be replaced with any other distance
    PU = I-UC*pinv(UC'*UC)*UC';
    maximum = 0;

    % Calculate the pixel most different from the already selected ones according to
    % the orthogonal projection (or any other distance selected)
    for n = 1:ns
        for m = 1:nl
            r = squeeze(HIM(n,m,:));
            result = PU*r;
            val = result'*result;
            if (val > maximum) 
                maximum = val;
                posx = n; posy = m;
            end
        end
    end
    
    % The next chosen pixel is the most different from the already chosen ones
    ti = squeeze(HIM(posx,posy,:));
    
    % Show positions of the above mentioned pixel at screen
    fprintf(1, '. found pixel @ coordinates x=%5d & y=%5d \n',posx,posy);

    % Store positions in the matrix of positions
    P(i+1,1)=posx; P(i+1,2)=posy;

    % Store positions in the matrix of pure pixels
    U = [U ti];

    % Visualize pixel selected on screen
    drawnow;
    text(posy,posx,'o','Margin',1,'HorizontalAlignment','center','FontSize',22,'FontWeight','light','FontName','Garamond','Color','yellow');
%     text(posy,posx,'o','Color','yellow');
end

% Obtener tiempo CPU actual | get current CPU time
t2 = cputime;

% Show total execution time of the algorithm
fprintf(1, '. Total CPU processing time .................... %6.3f [s] \n', (t2-t1) );
disp(' === Eng OSP ===');
