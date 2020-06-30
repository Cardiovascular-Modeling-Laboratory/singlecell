% step 1: translate P0 and P4 to the horizontal axis with the transformed
% P0 being at the origin
% translate the control points
P0_t = P0-P0;
P4_t = P4-P0;

% rotate the new control points 
%   determine the rotation angle: the angle P4_t forms with the positive
%   x-axis. we assume that theta>0 if {P4_t}_y>0 and theta<0 if {P4_t}_y<0.
%   then the rotation angle will be phi=-theta.
P4_t_x = P4_t(1);
P4_t_y = P4_t(2);
if P4_t_x==0 && P4_t_y>0 %if P4_t lies along positive y-axis
    theta=pi/2;
elseif P4_t_x==0 && P4_t_y<0 %if P4_t lies along negative y-axis
    theta=-pi/2;
elseif P4_t_y==0 %if P4_t lies along x-axis
    theta=0;
elseif P4_t_x>0 && P4_t_y>0 %if P4_t is strictly in quadrant 1
    theta=atan(abs(P4_t_y/P4_t_x));
elseif P4_t_x<0 && P4_t_y>0 %if P4_t is strictly in quadrant 2
    theta=pi-atan(abs(P4_t_y/P4_t_x));
elseif P4_t_x<0 && P4_t_y<0 %if P4_t is strictly in quadrant 3
    theta=atan(abs(P4_t_y/P4_t_x))-pi;
elseif P4_t_x>0 && P4_t_y<0 %if P4_t is strictly in quadrant 4
    theta=-atan(abs(P4_t_y/P4_t_x));
end
phi = -theta;

% transform the control points
rot_matrix = [cos(phi) -sin(phi); sin(phi) cos(phi) ];
P0_tt=rot_matrix *P0_t';
P4_tt=rot_matrix *P4_t';

% need to also transform the vectors. If we know the starting points and
% the vectors, then we know the ending point satisfies v=P_end-P_start,
% i.e, P_end=P_start+v. We can now apply the transformation to the endpoint
% and then determine the corresponding transformed vector using
% v_t=P_end_t-P_start_t
P0_end=P0+W1; % determine the initial endpoints
P4_end=P4+W2;

P0_end_t = P0_end-P0; % translate the endpoints
P4_end_t = P4_end-P0;
P0_end_tt=rot_matrix *P0_end_t'; % rotate the translated endpoints
P4_end_tt=rot_matrix *P4_end_t';

W1_tt=P0_end_tt-P0_tt; % use transformed endpoints to determine transformed vectors
W2_tt=P4_end_tt-P4_tt;

% to determine whether the vectors point towards each other: consider the
% midpoint of the transformed starting points. Translate the 2
% starting/ending points horizontally so that the 2 starting points lie on
% opposite sides of the vertical axis. we can then consider the lines
% connecting each starting point to its corresponding endpoint. These lines
% will intersect the vertical axis at either positive or negative values.
% If the connection is a desirable connection, then the y-intercepts of
% both lines will be positive. 


