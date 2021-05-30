function [ xy , K ] = InterfaceCurvature( I , v , NN , SR )

    Flag = 1;  % Flag = 1 --> Display plots
               % Flag = 0 --> Do not display plots

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Description of I/O variables:                                       %
    %                                                                     %
    % INPUT:-                                                             %
    %                                                                     %
    % I   : image matrix, can be 2D or 3D (RGB)                           %
    % v   : contour / level-set value, which should be in the range [0,1] %
    % NN  : number of neighboring points used to fit a circle             %
    % SR  : sub-sampling rate (should be integer > 0)                     %
    %                                                                     %
    % OUTPUT:-                                                            %
    %                                                                     %
    % K   : curvature vector (cell of dim equal to number of contours)    %
    % xy  : coordinates of the points on the contour                      %
    %                                                                     %
    %                                                                     %
    % NOTE: For reliable outputs, the image must be atleast 300x300       %
    %       and each contour must have atleast 500 points.                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Author: Jaya Kumar. A
    % E-mail: jkumar.res@gmail.com
    %
    % See the accompanying "doc.pdf" for algorithm details.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Formating the data %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Converting a 3D color image to 2D grayscale image
    
    if length(size(I))==3
        I = rgb2gray(I);
    end

    I = double(I);
    Imin = min(I(:));
    Imax = max(I(:));

    I = (I-Imin) / (Imax-Imin);  % Normalizing the data [0,1]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Padding the data to take care of   %%%
    %%% interfaces that hit the boundaries %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [Lx,Ly] = size(I);

    I2 = zeros(Lx+6,Ly+6);
    I2( 4:Lx+3 , 4:Ly+3 ) = I;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Contour extraction %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    N = NN;
    figure;
    subplot(2,2,1)
    xy = contour( I2 , [v v] );
    title('Contours');
    xlabel('X');
    ylabel('Y');
    axis equal;
    axis tight;
    %set(gca,'YDir','reverse');

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Conttour splitting %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Np = length( xy(1,:) );  
    Lc = xy(2,1);   % Number of points on the 1st contour
    ll = Lc+1;
    iLc = 1;
    
    while( ll < Np )
        Lc = [Lc ; xy(2,ll+1)];
        iLc = [iLc ; ll+1];
        ll = sum( Lc+1 );        
    end
    
    Nc = length( Lc );    % Number of contours
    
    
    % Cells to handle multiple contours with different sets of points
    
    x = cell(Nc,1);    % X-coordinates of the contour
    y = cell(Nc,1);    % Y-coordinates of the contour
    r = cell(Nc,1);    % Radius of curvature 
    rr = cell(Nc,1);   % 
    x2 = cell(Nc,1);   % Sampled X-coordinates
    y2 = cell(Nc,1);   % Sampled Y-coordinates
    a = cell(Nc,1);    % X-coordinate of the circle center
    b = cell(Nc,1);    % Y-coordinate of the circle center
    
    for i=1:Nc
        x{i} = xy( 1 , iLc(i)+1 : iLc(i)+Lc(i) )';
        y{i} = xy( 2 , iLc(i)+1 : iLc(i)+Lc(i) )';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Calculation of curvature %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=1:Nc
        %%%%%%%%%%%%%%%%%%%%
        %%% Sub-sampling %%%
        %%%%%%%%%%%%%%%%%%%%
        
        Np = length(x{i});

        ll = floor( Np/SR );
        xx = reshape( x{i}(1:SR*ll) , [SR ll] );
        yy = reshape( y{i}(1:SR*ll) , [SR ll] );

        xx = xx(1,:)';
        yy = yy(1,:)';
        Np = ll;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Locally fitting circles %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        r{i} = zeros(Np,1);
        a{i} = zeros(Np,1);
        b{i} = zeros(Np,1);

        for j=1:Np
            %%% Maintain contour periodicity
            j1 = mod( j-NN-1 , Np ) + 1;
            j2 = mod( j+NN-1 , Np ) + 1;
            
            if( j1 < j2 )
                [a{i}(j),b{i}(j),r{i}(j)] = FitCircle(xx( j1:j2 ),yy( j1:j2 ));
            else
                [a{i}(j),b{i}(j),r{i}(j)] = FitCircle( [xx( j1:end ) ; xx(1:j2)]...
                                          ,[yy( j1:end ) ; yy(1:j2)]);
            end
        end
                    
        rr{i} = abs(r{i});
        
        x2{i} = xx;
        y2{i} = yy;
        a{i} = real(a{i});
        b{i} = real(b{i});
    end
    
    if( Flag==1 )
        subplot(2,2,2);
        plot(rr{1});
        title('Radius of curvature');
        hold on;
        leg = cell(Nc,1);
        leg{1} = 'Contour1';

        for i=2:Nc
            plot(abs(rr{i}));
            leg{i} = sprintf('Contour%d',i);
        end
    
        ylim( [0 , sqrt(Lx^2+Ly^2)] );
        xlabel('Contour length');
        ylabel('R');
        legend( leg );
    
        subplot(2,2,4);
        plot3(x2{1},y2{1},rr{1},'-*');
        title('3D plot of radius of curvature');
        hold on;
    
        for i=2:Nc
            plot3(x2{i},y2{i},rr{i},'-*');
        end        
    
        box on;
        zlim( [0 , sqrt(Lx^2+Ly^2)] );
        xlabel('X');
        ylabel('Y');
        zlabel('R');
    end
    
    v = cell(Nc,1);
    K = cell(Nc,1);
    
    for i=1:Nc
        v{i} = [(x2{i}-a{i}) (y2{i}-b{i})];
        vv = sqrt( v{i}(:,1).^2 + v{i}(:,2).^2 );
        v{i} = v{i} ./ [vv vv];
        
        K{i} = v{i} ./ [rr{i} rr{i}];
    end
    
    if( Flag==1 )   
        subplot(2,2,3);
        [m,n] = size(I2);
        I3 = zeros(m,n,3);
        I3(:,:,1) = I2;
        I3(:,:,2) = I2;
        I3(:,:,3) = I2;
        %imagesc(I3);
        hold on;
            
        for i=1:Nc
            plot(x2{i},y2{i});
            h = quiver(x2{i},y2{i},K{i}(:,1),K{i}(:,2));
            set(h,'AutoScale','on', 'AutoScaleFactor', 2)
        end
    
        axis equal;
        axis tight;
        xlim([0 Ly]);
        ylim([0 Lx]);
        title('Negative curvature vector field');
        xlabel('X');
        ylabel('Y');
    end
    
    xy = cell(Nc,1);
    
    for i=1:Nc
        xy{i} = [x2{i} y2{i}];
    end
    
end