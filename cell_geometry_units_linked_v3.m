% File name: cell_geometry_units_linked_v3.m
%
% Last updated (by Will): 1/15/2019
% First version finished: 08/21/2009
%
% written by Anya Grosberg
% Disease Biophysics Group
% School of Engineering and Applied Sciences
% Harvard University, Cambridge, MA 02138
%
% The purpose of this file is twofold: First, we will be able to designate
% a cell shape and discretize the shape for simulation. To do this, we'll
% create a set of r vectors that describe the area of the cell as well as
% all the points within the cell for which the model will be run. Secondly,
% to compute the distance terms that will be needed during the simulation.
% This includes the general distance r-r' and the distance term within the
% biased diffusion equation (computed for every point r and every point r')
%
%
% Input:  1. Approximate total number of points (Num_points) from parameter file -- this number will be adjusted so that 
% there is a whole number of points used for the specified geometry
%         2. Area of the cell (A) from parameter file
%         - The user will be asked to choose the shape, and some properties
%         of the shape if approprate. Here are the current list of shapes
%         available:
%           [1] = Square
%           [2] = Rectangle; The user will then be asked to input an aspect
%           ratio
%           [3] = Circle
%           [4] = Triangle; The user will then be asked to input the angle
%           at the "top" of the triangle. It is assumed that the triangle
%           being considered is either EQUILATERAL or ISOSCELES so that only
%           one angle is needed to create the geometry
%           [5] = Stair, Webbed; (not adjusted from original code since it
%           was primarily used to fit parameters)
%           [6] = Stair, Plain; (not adjusted from original code since it
%           was primarily used to fit parameters)
%           [7] = Oval; The user will then be asked to input an aspect ratio
%           [8] = Concave shape: rectangle with triangular cutout; the user
%           will be asked to input a rectangular aspect ratio and an angle
%           for the triangle 
%           [9] = T shaped cell: a test case for building the network and
%           fitting some parameters. 
%           [10] = V shaped cell: another concave cell to test
%
% Global: 1. Fig_Num number of the previous figure number
%
%
% Output: 1. A set of vectors that describe the cell area space (mat_r)
%         2. Total number of points used - this is the same as length of r,
%         but it's easier to pass it on (Npts_t)
%         3. The differences in x coordinates - each row - vary r', column
%         - vary r (drx)
%         4. The differences in y coordinates - each row - vary r', column
%         - vary r (dry)
%         5. The unit area (dA)
%         6. The square of the distance between each r-r' pair
%         (dr_dist_squared)
%         7. The square of the distance between r and line r''-r'
%         (dist_to_line_sq)
%         8. Shape name - the name of the shape - this will be used for
%         saving the output files, so that when multiple shapes are run
%         from the same paramter file, the results don't overwrite
%         themselves (shape_name)
%         9. The distance between each pair |r'-r''| (dist_pair)
%        10. drx_norm -- normalized x distance between r and r'
%        11. dry_norm -- normalized y distance between r and r'
%        12. A list of indicies of mat_r that indicate which of the
%        elements should be zero. This list is equal to 0 for all convex
%        shapes, and contains the indicies of the coordinates of the points
%        outside the concave shapes (Concave_ind)
%        13. outline - is the coordinates of line segments necessary to
%        outline the shpae.

function [mat_r,Npts_t,drx,dry,dA,dr_dist_squared, dist_to_line_sq, shape_name,dist_pair,drx_norm,dry_norm,Concave_ind,outline,nuc_x,nuc_y,nuc_cx,nuc_cy,choice, outside_segs, inside_segs,bdry_mat,bdry_pts,out_ind,Num_points,A,nuc_radius]=cell_geometry_units_linked_v3(Num_points,A,nuc_radius)

%Declare a variable Fig_Num to be global - this represents the figure
%number, so that each set of figures generates in a different window.
global Fig_Num
Fig_Num=0;
%
%Ask the user which shape they would like to model
choice = input(...
    'What kind of shape would you like to model? 1==[square], 2==[rectangle], 3==[circle],\n 4==[triangle], 5==[stair, webbed], 6==[stair, plain],  7==[oval] ,\n 8==[concave rectangle]   9==[T shape],   10==[V shape]: ');

%Initilize the variable that will list all the concave shape empty indicies
%- this will be a vector the same length as mat_r or length 1 for
%non-concave shapes. The vector will have zeros in the indicies where there
%should be no integrin biding.
Concave_ind = 1;

%%
%Square geometry
if choice == 1
    l=sqrt(A); % length of the side of the square
    num_points_x = ceil(sqrt(Num_points)); %number of points in the x direction
    num_points_y = num_points_x; %fort the square these would be the same
    r_x = linspace(0,l,num_points_x)'; %generate the x coordinates (column vector)
    r_y = linspace(0,l,num_points_y); %generate the y coordinats (row vector)
    mat_r(:,1) = repmat(r_x,num_points_y,1);%write out x coordinates
    r_y_vec_temp = repmat(r_y,num_points_x,1); %generate the repeating column array of y vectors
    mat_r(:,2) = r_y_vec_temp(:); %set the y values into a single column
    shape_name = 'square'; %the name of the shape
    outline = [0 sqrt(A) sqrt(A) 0 0; 0 0 sqrt(A) sqrt(A) 0];
end
%%
%Rectangular geometry
if choice == 2
    %ask the user for the aspect ratio
    aspect_ratio = input('Please specify the aspect ratio, i.e. x/y: ');
    l=sqrt(A*aspect_ratio);
    w=l/aspect_ratio;
    num_points_x = ceil(sqrt(aspect_ratio*Num_points)); %number of points in the x direction
    num_points_y = ceil(num_points_x/aspect_ratio); %use aspect ratio to backtrack the number of points -- err on the side of more points
    r_x = linspace(0,l,num_points_x)'; %generate the x coordinates (column vector)
    r_y = linspace(0,w,num_points_y); %generate the y coordinats (row vector)
    mat_r(:,1) = repmat(r_x,num_points_y,1);%write out x coordinates
    r_y_vec_temp = repmat(r_y,num_points_x,1); %generate the repeating column array of y vectors
    mat_r(:,2) = r_y_vec_temp(:); %set the y values into a single column
    shape_name = ['rectangle_ar=' num2str(aspect_ratio)];
    outline = [0 sqrt(A*aspect_ratio) sqrt(A*aspect_ratio) 0 0; 0 0 sqrt(A/aspect_ratio) sqrt(A/aspect_ratio) 0]; %the coorinates for drawing the outline of the shape
end
%%
%Circular geometry - the best way to do this is to create the
%vectors in polar coordinates and then backtrack to calculate the
%cartesian coordinates, this way all the points will be spaced
%evenly - the number of points for each r will be proportional to r
if choice == 3
    r_max = sqrt(A/pi()); % radius that describes the circle
    num_points_r = ceil(sqrt(Num_points/pi)); %number of points in the r coordiante
    r_coord = linspace(0,r_max,num_points_r); % radii points
    num_points_t_max = num_points_r*2*pi(); %number of points in the outer circle - designed so that the density is the same as in the r direction
    n_theta=linspace(1,num_points_t_max,num_points_r);
    num_points_t = ceil(n_theta); %number of points in the theta coordinate for each r coordinate    
    dtheta_r = (2.*pi().*ones(1,num_points_r))./(num_points_t); %the step in theta for each r
    mat_r = [];
    for i=1:num_points_r %cycle through each r
        theta = linspace(0,2*pi()-dtheta_r(i),num_points_t(i)); %the distribution of theta
        r_x = r_coord(i).*cos(theta); %x coord
        r_y = r_coord(i).*sin(theta); %y coord
        mat_r = vertcat(mat_r,[r_x',r_y']); %populate the array with the r vectors
    end
    shape_name = 'circle';
    outline = [0 0 r_max]; %the coorinates for drawing the outline of the shape - in the circular case just the coordinates of the origin and outer radii
end
%%
%Triangular geometry, the user will be asked for the top angle
if choice == 4
    angle = input('Enter the angle at the top of this isosceles triangle in degrees: ');
    l = 2*sqrt(A)*sqrt(tand(angle/2)); %the length of the base triangle
    h = 2*A/l; %the height of the triangle
    num_points_base= ceil(sqrt(2*Num_points)); %the number of points in the base
    num_points_height = num_points_base; %the number of points in the height

    h_k = linspace(0,h,num_points_height); %y coordinate of each row, y=0 at the base
    length_hor_lines = l*(h-h_k)./h ; %the width of the triangle at each h_k
    num_points_hor = ceil((num_points_base/l).*length_hor_lines+1); %the number of points at each horizontal line;
    x_k1=l.*h_k./(2*h);
    mat_r =[]; %initilize the coordinate matrix
    for i=1:num_points_height %cycle through each height point starting from the bottom
        start_x = x_k1(i); %where the x coordinate starts
        r_x = linspace(start_x,start_x+length_hor_lines(i),num_points_hor(i))'; %the x coordinates
        mat_r = vertcat(mat_r,[r_x,repmat(h_k(i),num_points_hor(i),1)]); %populate the x,y coordinates
    end
    shape_name = ['triangle_' num2str(angle)];
    outline = [0 2*max(r_x) start_x 0; 0 0 h_k(num_points_height) 0]; %the coorinates for drawing the outline of the shape
end
%%
%Stair shape geometry; Includes webbing (i.e. there are two areas \Omega and \Omega'
if choice == 5
    Num_points=150;
    %At this point we're going to assume the dimensions of the stair shape
    %that correspond to the experiment, this can be edited at some other
    %point 
    A = 2083*10^(-12); % Area of the cell written here in square meters, |Delta| in the notes
    nuc_area = 129.4192*10^(-12); % The area of the nucleus in m 
    nuc_radius=sqrt(nuc_area/pi);
    
    x_length = 87.3e-6;% 85.5e-6;
    y_length = 37.9e-6; %35e-6;
    x_step1 = 28.5e-6; % 28e-6;
    y_step1 = 17.1e-6; % 16.5e-6;
    x_step2 = 58e-6; % 63e-6;
    y_step2 = 21.5e-6; % 12.5e-6;    
    aspect_ratio = x_length/y_length;
    num_points_x = ceil(sqrt(aspect_ratio*Num_points)); %number of points in the x direction
    num_points_y = ceil(num_points_x/aspect_ratio); %use aspect ratio to backtrack the number of points -- err on the side of more points

    r_x = linspace(0,x_length,num_points_x)'; %generate the x coordinates (column vector)
    r_y = linspace(0,y_length,num_points_y); %generate the y coordinats (row vector)
    mat_r(:,1) = repmat(r_x,num_points_y,1);%write out x coordinates
    r_y_vec_temp = repmat(r_y,num_points_x,1); %generate the repeating column array of y vectors
    mat_r(:,2) = r_y_vec_temp(:); %set the y values into a single column
    
    [~,idx2]=min(abs(r_x-x_step2));
    x_step2_approx=r_x(idx2);
    [~,idx1]=min(abs(r_x-x_step1));
    x_step1_approx=r_x(idx1);
    [~,idx2]=min(abs(r_y-y_step2));
    y_step2_approx=r_y(idx2);
    [~,idx1]=min(abs(r_y-y_step1));
    y_step1_approx=r_y(idx1);
    
    %Concave_ind
    Concave_ind = ones(size(mat_r(:,1)))-cast(((mat_r(:,1)<x_step1_approx & mat_r(:,2)>y_step1_approx) | (mat_r(:,1)>x_step2_approx & mat_r(:,2)<y_step2_approx)),'double');
    Concave_ind_act = find(Concave_ind);
    shape_name = 'stair_web';
    outline = [0 x_step2 x_step2 x_length x_length x_step1 x_step1 0 0; 0 0 y_step2 y_step2 y_length y_length y_step1 y_step1 0]; %the coorinates for drawing the outline of the shape
end
%%
%Stair shape geometry: no webbing
if choice == 6
    Num_points=150;
    %At this point we're going to assume the dimensions of the stair shape
    %that correspond to the experiment, this can be edited at some other
    %point. 
    A = 2083*10^(-12); % Approximate area of the cell written here in square meters, |Delta| in the notes
    nuc_area = 129.4192*10^(-12); % The area of the nucleus in m 
    nuc_radius=sqrt(nuc_area/pi);    
    
    x_length = 87.3e-6;% 85.5e-6;
    y_length = 37.9e-6; %35e-6;
    x_step1 = 28.5e-6; % 28e-6;
    y_step1 = 17.1e-6; % 16.5e-6;
    x_step2 = 58e-6; % 63e-6;
    y_step2 = 21.5e-6; % 12.5e-6;    
    aspect_ratio = x_length/y_length;
    
    %adjust the number of points so that the total number of points inside
    %the cell works with the user selection
    empty_area = ((x_length-x_step2)*y_step2 + x_step1*(y_length-y_step1));
    Num_points_adjust = Num_points*(1+empty_area/(x_length*y_length));
    num_points_x = ceil(sqrt(aspect_ratio*Num_points_adjust)); %number of points in the x direction
    num_points_y = ceil(num_points_x/aspect_ratio); %use aspect ratio to backtrack the number of points -- err on the side of more points

    r_x = linspace(0,x_length,num_points_x)'; %generate the x coordinates (column vector)
    r_y = linspace(0,y_length,num_points_y); %generate the y coordinats (row vector)
    mat_r_temp(:,1) = repmat(r_x,num_points_y,1);%write out x coordinates
    r_y_vec_temp = repmat(r_y,num_points_x,1); %generate the repeating column array of y vectors
    mat_r_temp(:,2) = r_y_vec_temp(:); %set the y values into a single column
    
    [~,idx2]=min(abs(r_x-x_step2));
    x_step2_approx=r_x(idx2);
    [~,idx1]=min(abs(r_x-x_step1));
    x_step1_approx=r_x(idx1);
    [~,idx2]=min(abs(r_y-y_step2));
    y_step2_approx=r_y(idx2);
    [~,idx1]=min(abs(r_y-y_step1));
    y_step1_approx=r_y(idx1);
    
    cell_ind = (mat_r_temp(:,1)<=x_step2_approx & mat_r_temp(:,2)<=y_step1_approx) | ((mat_r_temp(:,1)>=x_step1_approx & mat_r_temp(:,1)<=x_step2_approx) & ...
        (mat_r_temp(:,2)>=y_step1_approx & mat_r_temp(:,2) <= y_step2_approx)) | (mat_r_temp(:,1) >= x_step1_approx & mat_r_temp(:,2) >= y_step2_approx);
    mat_r(:,1) = mat_r_temp(find(cell_ind),1);
    mat_r(:,2) = mat_r_temp(find(cell_ind),2);

    cell_ind_1 = mat_r(:,2)<y_step2_approx;
    cell_ind_2 = mat_r(:,1)>x_step2_approx;
    x_step2 = max(mat_r(find(cell_ind_1),1));
    y_step2 = min(mat_r(find(cell_ind_2),2));
    cell_ind_3 = mat_r(:,2)>y_step1_approx;
    cell_ind_4 = mat_r(:,1)<x_step1_approx;
    x_step1 = min(mat_r(find(cell_ind_3),1));
    y_step1 = max(mat_r(find(cell_ind_4),2));
    clear mat_r_temp cell_ind*

    shape_name = 'stair_plain';
    %the coorinates for drawing the outline of the shape
    outline = [0 x_step2 x_step2 x_length x_length x_step1 x_step1 0 0; 0 0 y_step2 y_step2 y_length y_length y_step1 y_step1 0];
end
%%
% Oval geometry
if choice==7
    %ask the user for the aspect ratio
    AR = input('Please specify the aspect ratio, i.e. x/y: ');    
    Area=A;
    % Determine a and b:
    % AR=a/b and Area=pi*a*b. So Area/b^2=pi*(a/b)=pi*AR, i.e.,
    % b=sqrt(Area/pi*AR). With b, we obtain a: Ar=a/b, i.e, a=b*AR
    b=sqrt(Area/(pi*AR));
    a=b*AR;

    % Check that a>b
    if a<=b
        disp(['a = ', num2str(a), ' and b = ', num2str(b), '. Need a>b. Input different values.'])
        return
    end

    disp(['Checkpoint 1 reached: a = ', num2str(a), ', b = ', num2str(b)]);

    % k=number of concentric ellipses to create
    Ne = ceil(sqrt(Num_points/pi));
    
    % define delta_r so that a-b is a constant multiple of delta_r if AR>=2
    % and a constant multiple of b is AR<2
    if AR<2    
        delta_r=b/Ne;
    else delta_r=(a-b)/Ne;
    end

    a_values = 0:delta_r:a; % vector of a values
    a_values=a-a_values;
    b_values = 0:delta_r:b; % vector of b values
    b_values=b-b_values;

    length_b=length(b_values);
    if b_values(end)==0
        b_values=b_values(1:length_b-1);
    %    A2=a_values(length(b_values)+1:end);
        a_values=a_values(1:length(b_values));
    else
    %    A2=a_values(length(b_values)+1:end);
        a_values=a_values(1:length(b_values));
    end

    % Right now, the elements of a_values and b_values are in decreasing order.
    % For the computation below, we want the elements to be in increasing
    % order. So flip them:
    a_k=fliplr(a_values);
    b_k=fliplr(b_values);

    ecc=sqrt(1-(b_k./a_k).^2); %vector of eccentricity for the concentric ellipses

    % If the eccentricity is too close to 1, placing points on the ellipse as
    % well as on the x-axis will create a high density of points near the
    % x-axis. To avoid this, set a tolerance and remove any ellipse with
    % eccentricity larger than this tolerance
    tol=input('Choose an eccentricity tolerance to rule out degenerate ellipse, i.e. 0.99 or 0.999, etc.: ');
    index=find(ecc>tol);
    a_k(index)=[];
    b_k(index)=[];
    ecc(index)=[];
    clear index

    % vector of arclenths for the concentric ellipses; matlab needs the 2nd input to be squared
    TotalArcLength=a_k.*ellipticE(2*pi,ecc.^2); 
    disp(['Checkpoint 2 reached: Perimeter of ellipse is ', num2str(TotalArcLength(end))]);

    %
    % max number of points on each ellipse
    K=linspace(Ne,2*Ne-1,size(b_k,2));
    n_k=K.*ellipticE(2*pi,ecc.^2);
    N_k=ceil(n_k);

    % Check that every element of N is even. If not, add 1 to make it an even
    % number
    if any(mod(N_k,2)==1) % checks if any element is an odd number
        index=find(mod(N_k,2));
        N_k(index)=N_k(index)+1;
    end
    clear index

    disp(['Checkpoint 3 reached: There will be ', num2str(length(N_k)), ' concentric ellipses']);

    % There are k non-degenerate concentric ellipses. For each ellipse, determine the
    % different parametric paramter t values (which corresponds to the number of points on each
    % ellipse)
    syms t
    T=zeros(sum(N_k),1); %initialize vector containing all parametric_t parameter values
    l=0;
    for j=1:length(N_k)
        param_t=zeros(1,N_k(j)); % initialize param_t vector
        for i=0:N_k(j)-1
            param_t_value=vpasolve(a_k(j)*ellipticE(t,ecc(j)^2)==i*TotalArcLength(j)/N_k(j),t);
            param_t(i+1)=param_t_value;
        end
        T(l+1:l+N_k(j),1) = param_t';
        l=N_k(j)+l;

        if mod(j,3)==0
            disp(['Number of points for ellipses 1 thru ', num2str(j),' determined']);
        end

    end
    T=double(T); % make T have "double" class instead of "sym" class

    disp('Checkpoint 4 reached: Calculating (x,y) coordinates');

    % Since points should be the same size as T, make A, B the
    % same size as T
    repeated_A=zeros(sum(N_k),1);
    repeated_B=zeros(sum(N_k),1);
    l=0;
    for j=1:length(N_k)
        repeated_A(l+1:l+N_k(j),1)=repmat(a_k(j),N_k(j),1);
        repeated_B(l+1:l+N_k(j),1)=repmat(b_k(j),N_k(j),1);
        l=N_k(j)+l;
    end
    a_k=repeated_A;
    b_k=repeated_B;

    % Now create points along the concentric ellipses
    x_points=a_k.*sin(T);
    y_points=b_k.*cos(T);

    % get positive axis points by looking at inner-most ellipse
    x_axis=x_points(1:N_k(1))';
    associated_points=y_points(1:N_k(1))';

    index=find(x_axis>0 & associated_points>0); % find index of points in 1st quadrant
    smallest_y=min(associated_points(index)); % of these points, find the smallest y value
    index2= associated_points(index)==smallest_y; % determine the index where smallest y value occurs
    index(index2)=[]; % remove this index from the "index" list

    pos_x_axis=x_points(index);
    pos_y_axis=zeros(1,length(pos_x_axis))';

    x_points=vertcat(x_points,pos_x_axis);
    y_points=vertcat(y_points,pos_y_axis);

    % insert negative axis points
    neg_x_axis=-pos_x_axis;
    x_points=vertcat(x_points,neg_x_axis);
    y_points=vertcat(y_points,pos_y_axis);

    % insert origin
    x_points=vertcat(x_points,0);
    y_points=vertcat(y_points,0);

    % Create matrix of (x,y) values
    mat_r=[x_points, y_points];

    shape_name = 'oval';
    temp_t=linspace(0,2*pi);
    X_ellipse_outline=a*cos(temp_t);
    Y_ellipse_outline=b*sin(temp_t);
    outline = [X_ellipse_outline; Y_ellipse_outline]; %the coorinates for drawing the outline of the shape

    disp('If you do not like the lattice point distribution used, consider either increasing tolerance or decreasing Num_points')
end
%%
% Concave geometry: square with trianglular cutout
if choice==8
    aspect_ratio= input('Please specify the aspect ratio, i.e. x/y: ');    
    angle= input('Enter the angle at the top of the embedded isosceles triangle (in degrees): ');
   
    checkpt=1-aspect_ratio/(4*tand(angle/2));
    if checkpt<=0
        disp('This choice of aspect ratio and triangle angle are incompatible with the shape requirements. Pick different values.')
        return
    end
    
    A_rect=A/checkpt;
    A_tri=A_rect-A;
    
    l=sqrt(A_rect*aspect_ratio);
    w=l/aspect_ratio;
    h=2*A_tri/l;
    
    % the discretization is based on previous discretization techniques.
    % the triangular portion has it's peak placed at the origin. Due to the
    % symmetry of the shape, we'll only construct 2 regions (one rectangle,
    % one triangle) and then reflect them about the y-axis to generate the
    % other regions of the shape
    L=l/2;
    W=w-h;
    num_points_x = ceil(sqrt(aspect_ratio*Num_points/4)); %number of points in the x direction
    num_points_y = ceil(num_points_x/aspect_ratio); %use aspect ratio to backtrack the number of points -- err on the side of more points
    r_x = linspace(0,L,num_points_x)'; %generate the x coordinates (column vector)
    r_y = linspace(0,W,num_points_y); %generate the y coordinats (row vector)
    mat_r_R1(:,1) = repmat(r_x,num_points_y,1);%write out x coordinates
    r_y_vec_temp = repmat(r_y,num_points_x,1); %generate the repeating column array of y vectors
    mat_r_R1(:,2) = r_y_vec_temp(:); %set the y values into a single column
    
    mat_r_R1p=mat_r_R1;
    mat_r_R1p(:,1)=-mat_r_R1p(:,1);
    
    num_points_base= ceil(sqrt(2*Num_points/4)); %the number of points in the base
    num_points_height = ceil(num_points_base/2); %num_points_base; %the number of points in the height
    h_k = linspace(0,-h,num_points_height); %y coordinate of each row, y=0 at the base
    length_hor_lines = L.*(1+h_k./h) ; %the width of the triangle at each h_k
    num_points_hor = ceil((num_points_base/L).*length_hor_lines+1); %the number of points at each horizontal line;
    x_k1=-L.*h_k./h;
    mat_r_R2 =[]; %initilize the coordinate matrix
    for i=2:num_points_height %cycle through each height point starting from just above the triangle base
        %start_x = h_k(i)/slope; %where the x coordinate starts
        start_x = x_k1(i); %where the x coordinate starts
        r_x = linspace(start_x,start_x+length_hor_lines(i),num_points_hor(i))'; %the x coordinates
        mat_r_R2 = vertcat(mat_r_R2,[r_x,repmat(h_k(i),num_points_hor(i),1)]); %populate the x,y coordinates
    end
    
    mat_r_R2p=mat_r_R2;
    mat_r_R2p(:,1)=-mat_r_R2p(:,1);
    
    % combine the lattices for the different regions into one region and
    % remove duplicate points
    mat_r=unique([mat_r_R1; mat_r_R1p; mat_r_R2; mat_r_R2p],'rows');
    
    % identify the points that outline the shape
    p1=[min(mat_r_R2p(:,1)), min(mat_r_R2p(:,2))];
    p2=[0, 0];
    p3=[max(mat_r_R2(:,1)), min(mat_r_R2(:,2))];
    p4=[max(mat_r_R1(:,1)), max(mat_r_R1(:,2))];
    p5=[min(mat_r_R1p(:,1)), max(mat_r_R1p(:,2))];
    outline_temp=[p1;p2;p3;p4;p5;p1];
    outline=outline_temp';
    
    shape_name = ['concave_ar=' num2str(aspect_ratio) '_theta=' num2str(angle)];
    
end
%%
% Concave geometry: T shaped cell (this shape has fixed area and nucleus size
% since it's primarily for general purpose testing and not for validation).
% Shape dimensions are taken from reference paper Pathak2007
if choice == 9
%     A=480e-12; % ~cell area 
%     nuc_radius=5e-6;
%     
%     L1=6; % cell length1, in micrometers 
%     H=40; % cell height, in micrometers 
%     L2=46; % cell length2, in micrometers 

    L1=12; % cell length1, in micrometers 
    H=(A*10^12+L1^2)/((46/40+1)*L1); % cell height, in micrometers 
    L2=(46/40)*H; % cell length2, in micrometers 

    l=L1/2;
    w=H;
    aspect_ratio=l/w;
    num_points_x = 2; %number of points in the x direction
    num_points_y = ceil(num_points_x/aspect_ratio); %use aspect ratio to backtrack the number of points -- err on the side of more points
    
    r_x = linspace(0,l,num_points_x)'; %generate the x coordinates (column vector)
    r_y = linspace(0,w,num_points_y); %generate the y coordinats (row vector)
    mat_r_R1(:,1) = repmat(r_x,num_points_y,1);%write out x coordinates
    r_y_vec_temp = repmat(r_y,num_points_x,1); %generate the repeating column array of y vectors
    mat_r_R1(:,2) = r_y_vec_temp(:); %set the y values into a single column
    
    clear l w aspect_ratio num_points_x num_points_y
    dy=r_y(2)-r_y(1);
    dx=r_x(2)-r_x(1);
    
    r_y2=r_y(min(find(r_y>=H-L1-1))):dy:max(r_y);
    r_x2=0:dx:L2/2;
    r_x2=r_x2';
    num_points_x = numel(r_x2); %number of points in the x direction
    num_points_y = numel(r_y2); %use aspect ratio to backtrack the number of points -- err on the side of more points
    mat_r_R2(:,1) = repmat(r_x2,num_points_y,1);%write out x coordinates
    r_y_vec_temp = repmat(r_y2,num_points_x,1); %generate the repeating column array of y vectors
    mat_r_R2(:,2) = r_y_vec_temp(:); %set the y values into a single column

    mat_r_R1p(:,1)=-mat_r_R1(:,1);
    mat_r_R1p(:,2)=mat_r_R1(:,2);
    mat_r_R2p(:,1)=-mat_r_R2(:,1);
    mat_r_R2p(:,2)=mat_r_R2(:,2);
    
    mat_r=[mat_r_R1; mat_r_R1p; mat_r_R2; mat_r_R2p];
    mat_r=unique(mat_r,'rows');
    xs=mat_r(find(mat_r(:,2)==0),1);

    outline_x=[max(xs) max(xs) max(r_x2) max(r_x2) -max(r_x2) -max(r_x2) -max(xs) -max(xs) max(xs)];
    outline_y=[min(r_y) min(r_y2) min(r_y2) max(r_y2) max(r_y2) min(r_y2) min(r_y2) min(r_y) min(r_y)];
    outline=[outline_x; outline_y];

    mat_r=mat_r.*10^(-6);
    outline=outline.*10^(-6);
    shape_name = 'Tshape';
end
%%
% Concave geometry: V shaped cell (this shape has fixed area and nucleus size
% since it's primarily for general purpose testing and not for validation).
% Shape dimensions are taken from reference paper Pathak2007
if choice == 10
    A=468e-12; % ~cell area
    d=6; % width of portion of cell
    N_L=4; % number of lines in half the V shape that will be discretized
    ymax=sqrt(46^2-23^2);
    m=sqrt(46^2-23^2)/23;
    mat_r_R1=[];
    xs=[];
    ys=[];
    
    for i=1:N_L
        d_i=(i-1)*d/(N_L-1);
        b_i=(23/46 + (46^2-23^2)/(23*46))*d_i;
        xi_max=(ymax-b_i)/m;
        
        num_points_base= ceil(sqrt((N_L-i+1)*(Num_points)));
        
        Xi=linspace(0,xi_max,num_points_base);
        Yi=m.*Xi+b_i;
        
        mat_r_R1=[mat_r_R1; [Xi' Yi']];
        
        if i==1 || i==N_L
            xs= [xs, min(Xi), max(Xi)];
            ys= [ys, min(Yi), max(Yi)];
        end
    end

    mat_r_R2(:,1)=-mat_r_R1(:,1);
    mat_r_R2(:,2)=mat_r_R1(:,2);
    
    mat_r=[mat_r_R1; mat_r_R2];
    mat_r=unique(mat_r,'rows');


    outline_x=[xs(1) xs(2) xs(4) xs(3) -xs(4) -xs(2) xs(1)];
    outline_y=[ys(1) ys(2) ys(4) ys(3) ys(4) ys(2) ys(1)];
    outline=[outline_x; outline_y];

    mat_r=mat_r.*10^(-6);
    outline=outline.*10^(-6);
    shape_name = 'Vshape';
end
%% 
if(choice<1 || choice>10 )
    %choice of non-exiting geometry
    disp('No such geomtry exists, stop the code (Ctrl^C) and re-run');
    return
end
%%
% Ask the user to place the nucleus in the cell
nuc_decide=input('Where do you want the nucleus? "Center"=1, "Randomly"=2, "Manually"=3: ');
if nuc_decide==1
    nuc_cx = (min(mat_r(:,1))+max(mat_r(:,1)))/2; %x-coordinate of center of nucleus
    nuc_cy = (min(mat_r(:,2))+max(mat_r(:,2)))/2; %y-coordinate of center of nucleus
elseif nuc_decide==2
    [nuc_cx, nuc_cy] = SymObjCtrProb_v3(A,choice,nuc_radius,mat_r);
else
    nuc_cx_temp = input(['Enter the x-coordinate for the center of the nucleus relative to the cell center. Assume the scale is micrometers.',...
        '\nFor example, "1.3" means 1.3um to the right of the cell center and "-5" means 5um to the left of the cell center: '] );
    nuc_cy_temp = input(['Enter the y-coordinate for the center of the nucleus relative to the cell center. Assume the scale is micrometers.',...
        '\nFor example, "1.3" means 1.3um above the cell center and "-5" means 5um below the cell center: '] );
    nuc_cx=nuc_cx_temp*10^(-6)+(min(mat_r(:,1))+max(mat_r(:,1)))/2;
    nuc_cy=nuc_cy_temp*10^(-6)+(min(mat_r(:,2))+max(mat_r(:,2)))/2;
end

clear t
t=linspace(0,1);
% nucleus curve
nuc_x=nuc_radius.*cos(2.*pi.*t)+nuc_cx;
nuc_y=nuc_radius.*sin(2.*pi.*t)+nuc_cy;
%%
%Calculate the total number of points that was used to make the r vectors
Npts_t = length(mat_r);
%Calculate the unit area
%The unit area needs to be adjusted by the number of points in the ECM
%island
if length(Concave_ind)~=1
    dA = A/(sum(Concave_ind));
else
    dA = A/Npts_t; %the unit area - total area=1
end
%%
%Calculate the r-r' vector for each pair - makes an array
r_col_x = repmat(mat_r(:,1),1,Npts_t); %make the array with each column having rx coorinates
r_row_x = r_col_x'; %make the array with each row having rx coorinates
r_col_y = repmat(mat_r(:,2),1,Npts_t); %make the array with each column having ry coorinates
r_row_y = r_col_y'; %make the array with each row having ry coorinates

if choice~=3 && choice~=7
    % relevent cell boundary geometry
    geo=outline'; 

    % identify all points that lie on the boundary of the cell. Also, determine 
    % which segments lie on the boundary. a 1 in the (i,j) position of bdry_mat
    % means the segment connecting point i and point j lies on the boundary of
    % the cell
    bdry_mat=zeros(size(r_col_x));
    t=linspace(0,1);
    mtx=unique(mat_r(:,1));
    mty=unique(mat_r(:,2));
    xd=min(abs(mtx(1:end-1,1)-mtx(2:end,1)))/2;
    yd=min(abs(mty(1:end-1,1)-mty(2:end,1)))/2;
    bdry_pts=[];
    for j=1:size(geo,1)-1
        P1=geo(j,:);
        P2=geo(j+1,:);
        seg_x=P1(1).*t+P2(1).*(1-t);
        seg_y=P1(2).*t+P2(2).*(1-t);
        Seg_x1=[seg_x'-xd seg_y'];
        Seg_x2=[seg_x'+xd seg_y'];
        Seg_y1=[seg_x' seg_y'-yd];
        Seg_y2=[seg_x' seg_y'+yd];
        Box=[Seg_x1;Seg_x2; Seg_y1; Seg_y2];
        k=convhull(Box(:,1),Box(:,2));
        Region=Box(k,:);
        in=inpoly(mat_r,Region);
        bdry_pt_idx=find(in==1);
        bdry_pts=[bdry_pts; mat_r(bdry_pt_idx,:)];

        meh=[nchoosek(bdry_pt_idx,2);fliplr(nchoosek(bdry_pt_idx,2))];
        bdry_mat(sub2ind(size(bdry_mat),meh(:,1),meh(:,2)))=1;
    end
    [out_ind,bdry_pts]=boundary_determ_v2(outline, mat_r,choice,bdry_pts);
else
    [out_ind,bdry_pts]=boundary_determ_v2(outline, mat_r,choice,[]);
    % relevent cell boundary geometry
    geo=bdry_pts; 

    % identify all points that lie on the boundary of the cell. Also, determine 
    % which segments lie on the boundary. a 1 in the (i,j) position of bdry_mat
    % means the segment connecting point i and point j lies on the boundary of
    % the cell
    bdry_mat=zeros(size(r_col_x));
    t=linspace(0,1);
    mtx=unique(mat_r(:,1));
    mty=unique(mat_r(:,2));
    xd=min(abs(mtx(1:end-1,1)-mtx(2:end,1)))/2;
    yd=min(abs(mty(1:end-1,1)-mty(2:end,1)))/2;
    for j=1:size(geo,1)-1
        P1=geo(j,:);
        P2=geo(j+1,:);
        seg_x=P1(1).*t+P2(1).*(1-t);
        seg_y=P1(2).*t+P2(2).*(1-t);
        Seg_x1=[seg_x'-xd seg_y'];
        Seg_x2=[seg_x'+xd seg_y'];
        Seg_y1=[seg_x' seg_y'-yd];
        Seg_y2=[seg_x' seg_y'+yd];
        Box=[Seg_x1;Seg_x2; Seg_y1; Seg_y2];
        k=convhull(Box(:,1),Box(:,2));
        Region=Box(k,:);
        in=inpoly(mat_r,Region);
        bdry_pt_idx=find(in==1);

        meh=[nchoosek(bdry_pt_idx,2);fliplr(nchoosek(bdry_pt_idx,2))];
        bdry_mat(sub2ind(size(bdry_mat),meh(:,1),meh(:,2)))=1;
    end    
end

% organize the boundary points so that you have a polygon that can be drawn
% as a closed contour
points.XY = bdry_pts;
P = tspo_ga(points);
rout=P.optRoute;
bdry_pts=[bdry_pts(rout,1), bdry_pts(rout,2)];

if choice==5 || choice==6 || choice==8 || choice==9 || choice==10 % if geometry is concave
    % identify region surrounding boundary
    bdry_reg_left=bdry_pts(:,1)-xd/2;
    bdry_reg_right=bdry_pts(:,1)+xd/2;
    bdry_reg_top=bdry_pts(:,2)+yd/2;
    bdry_reg_bottom=bdry_pts(:,2)-yd/2;
    br1=[bdry_reg_left, bdry_reg_bottom];
    br2=[bdry_reg_left, bdry_reg_top];
    br3=[bdry_reg_right, bdry_reg_bottom];
    br4=[bdry_reg_right, bdry_reg_top];
    bdry_reg1=[br1; br2; br3; br4];
    k=boundary(bdry_reg1(:,1),bdry_reg1(:,2),0.9);
    bdry_Region=bdry_reg1(k,:);
    
    if choice==9
        k=boundary(bdry_Region(:,1),bdry_Region(:,2),0.97);
        bdry_region2=[bdry_Region(k,1), bdry_Region(k,2)]; % the outter portion of the boundary region
    elseif choice==6 || choice==10
        bdry_region2=bdry_Region; % the outter portion of the boundary region
    end
    
    drx_temp = r_row_x - r_col_x; %the differences in x coordinates - each row - vary r', column - vary r
    dry_temp = r_row_y - r_col_y; %the differences in y coordinates - each row - vary r', column - vary r
    
    % this matrix indicates which segments do NOT lie on the boundary. so they
    % either lie completely on the interior of the cell or they will cross the
    % boundary
    inside_segs=zeros(size(bdry_mat));
    outside_segs=zeros(size(bdry_mat));
    seg_mat=1-bdry_mat; 
    t=linspace(0,1);
    
    for i1=1:size(seg_mat,1)
        pts=find(seg_mat(i1,:)>0);
        for i2=1:numel(pts)
            P1=mat_r(i1,:);
            P2=mat_r(pts(i2),:);

            % create segment connecting P1 to P2
            test_pts_x=P1(1).*t+P2(1).*(1-t);
            test_pts_y=P1(2).*t+P2(2).*(1-t);
            test_pts=[test_pts_x' test_pts_y'];
            test_pts(1,:)=[];
            test_pts(end,:)=[];

            % determine if segment is contained within the boundary geometry
            in=inpoly(test_pts,bdry_pts);
            in1=inpoly(test_pts,bdry_region2);

            % if segment is completely contained within the geometry or
            % completely outside the geometry, seg_mat(P1,P2)=1 and
            % seg_mat(P2,P1)=1. otherwise, seg_mat(P1,P2)=0 and
            % seg_mat(P2,P1)=0
            if sum(in==0)==98 % completely contained outside the geometry
                seg_mat(i1,pts(i2))=1;
                seg_mat(pts(i2),i1)=1;
                outside_segs(i1,pts(i2))=1;
                outside_segs(pts(i2),i1)=1;
            elseif sum(in==1)==98 % segment is completely contained within the geometry
                seg_mat(i1,pts(i2))=1;
                seg_mat(pts(i2),i1)=1;
                inside_segs(i1,pts(i2))=1;
                inside_segs(pts(i2),i1)=1;
            end
            if sum(in1==1)==98 % segment is inside geometry (redundancy check)
                seg_mat(i1,pts(i2))=1;
                seg_mat(pts(i2),i1)=1;
                inside_segs(i1,pts(i2))=1;
                inside_segs(pts(i2),i1)=1;                
            end
        end
    end    
    
    % create the matrix of relevant connections. (i,j)=1 if the segment
    % joining i and j is contained within the boundary of the cell or lies
    % completely outside the boundary of the cell. (i,j)=0 if not
    rel_connect=sign(inside_segs+outside_segs+bdry_mat); 

    drx = drx_temp.*rel_connect;
    dry = dry_temp.*rel_connect;
    
    clear x2_rx x1_rx y2_ry y1_ry fact2_x fact2_y fact1_x fact1_y Lines_cross_empty1 Lines_cross_empty2    
else
    drx = r_row_x - r_col_x; %the differences in x coordinates - each row - vary r', column - vary r
    dry = r_row_y - r_col_y; %the differences in y coordinates - each row - vary r', column - vary r
    inside_segs=ones(size(drx));
    outside_segs=zeros(size(drx));    
end

clear r_col_x r_col_y r_row_x r_row_y

%This means that to sum up over all r', I will need the command sum(drx,2)
%the square of the distance between each r-r' pair (note that this is a
%symmetric square matrix
dr_dist_squared = drx.^2+dry.^2; %the distance between each r'-r pair squared
dist_pair = sqrt(dr_dist_squared); %the distance between each pair of points

%%
%Calculate the distance between r and each line r''- r'
drx_trip_constr3 = repmat(drx,[1,1,Npts_t]); %repeat the differences in x into the third dimension that will be r'' (currently constant since only two variables), i.e. (r-r')x
dry_trip_constr3 = repmat(dry,[1,1,Npts_t]); %repeat the differences in y into the third dimension that will be r'' (currently constant since only two variables), i.e. (r-r')y
drx_trip_constr = permute(drx_trip_constr3,[3 2 1]); %flip so that the dimension along r is constant, i.e. (r''-r')x
dry_trip_constr = permute(dry_trip_constr3,[3 2 1]); %flip so that the dimension along r is constant, i.e. (r''-r')y
dr_dist_squaredr = permute(repmat(dr_dist_squared,[1,1,Npts_t]),[3 2 1]); % 3D array for |r''-r'|^2
dot_product_trip = drx_trip_constr.*drx_trip_constr3 + dry_trip_constr.*dry_trip_constr3;
to_line_x = drx_trip_constr3.*dr_dist_squaredr - drx_trip_constr.*dot_product_trip; % (r-r')|r''-r'|^2 - (r''-r')(dot product)
to_line_y = dry_trip_constr3.*dr_dist_squaredr - dry_trip_constr.*dot_product_trip;
dist_to_line_sq_temp = (to_line_x.^2 + to_line_y.^2)./(dr_dist_squaredr.^2); %the square of the distance from r to line (r''-r')
clear to_line_x to_line_y drx_trip_constr3 dry_trip_constr3 drx_trip_constr dry_trip_constr dot_product_trip
%the following is to get rid of the NaN from the distance - this appears
%when r''=r', and physically it should be zero at these point
temp_zero = zeros(size(dist_to_line_sq_temp));
dist_to_line_sq = dist_to_line_sq_temp;
dist_to_line_sq(find(isnan(dist_to_line_sq_temp))) = temp_zero(find(isnan(dist_to_line_sq_temp)));
clear dist_to_line_sq_temp temp_zero
%normalized drx and dry and get rid of NaN
temp_zero = zeros(size(drx));
drx_norm = drx./dist_pair;
drx_norm(find(isnan(drx_norm))) = temp_zero(find(isnan(drx_norm)));
temp_zero = zeros(size(dry));
dry_norm = dry./dist_pair;
dry_norm(find(isnan(dry_norm))) = temp_zero(find(isnan(dry_norm)));

Fig_Num = Fig_Num +1;
figure(Fig_Num)

hold on
%
tmp = size(outline);
hold on;
if tmp(1) ==1 %if circle; circle behaves weird with "patch" so make nucleus a filled circle
    rectangle('Position',[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
    rectangle('Position',[nuc_cx-nuc_radius,nuc_cy-nuc_radius,2*nuc_radius,2*nuc_radius],'Curvature',[1,1],'FaceColor','b');% this line plots the nucleus as a filled blue circle
else
    line(outline(1,:),outline(2,:),'Color',[0 0 0],'LineWidth',2);
    patch(nuc_x,nuc_y,'blue','facealpha',0.3); % this line plots the nucleus as a transparant blue circle
end
plot(mat_r(:,1),mat_r(:,2),'*r') %plot the shape for the user to check
if length(Concave_ind) ~= 1 %if stair shape with web then plot the empty points
    plot(mat_r(Concave_ind_act,1),mat_r(Concave_ind_act,2),'+b');
end
hold off;
axis equal;
set(gca,'Units','centimeters','OuterPosition',[5,5,6,6],'Position',[3,1,9.5,9.5]);
hold off
end