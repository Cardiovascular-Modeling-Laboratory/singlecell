% this code determines the tension for links over the concave region in a
% different way: The total force for the link i,j is defined as
% F_ij=(T*L0+EA*u_ij)e_ij where e_ij is the normalized vector pointing from
% i to j and u_ij=(d-L0)/L0. L0 is the initial length of the link, d is the
% distance between the points and T is the tension parameter
function [L_0_rel, d_rel,EA] = bdry_tension2(dist_pair,outside_segs,mat_r,outline)
% necessary parameters
EA=2000e-9; 
sig=0.7e-3;
l_f=EA/sig;

d_rel=dist_pair.*outside_segs;
R_rel=zeros(size(d_rel));
L_0_rel=zeros(size(d_rel));

for i1=1:size(d_rel,1)
    pts=find(d_rel(i1,:)>0);
    for i2=1:numel(pts)
        P0=mat_r(i1,:);
        P1=mat_r(pts(i2),:);
        d=d_rel(i1,pts(i2));

        r=(l_f/24)^(1/3)*d^(2/3);
        
        % determine circle center
        syms x0 y0
        S=vpasolve([ (P0(1)-x0)^2+(P0(2)-y0)^2==r^2, (P1(1)-x0)^2+(P1(2)-y0)^2==r^2],[x0 y0]);
        X0=S.x0;
        Y0=S.y0;

        % 2 possible centers: determine which center is correct (circle center can
        % be determined as follows: create a triangle with vertices at P0, P1, and
        % the center of the circle. The correct circle will have a center where the
        % trianglular shape only intersects with the boundary at points P0 and P1
        % while the incorrect circle will yield a triangluar shape that interesects
        % the boundary at more points than P0 and P1)
        xtest1=double(X0(1)); ytest1=double(Y0(1));
        pgon1=[P0;P1; xtest1 ytest1;P0];
        g1=unique(round(intersections(pgon1(:,1),pgon1(:,2),outline(1,:),outline(2,:)),10));
        xtest2=double(X0(2)); ytest2=double(Y0(2));
        pgon2=[P0;P1; xtest2 ytest2;P0];
        g2=unique(round(intersections(pgon2(:,1),pgon2(:,2),outline(1,:),outline(2,:)),10));

        if numel(g1)<=numel(g2)
            x0=double(X0(1));
            y0=double(Y0(1));
        else
            x0=double(X0(2));
            y0=double(Y0(2));
        end

        % determine arclength of the circular arc using s=R*theta where
        % cos(theta)=dot(A,B)/norm(A)norm(B).
        A=P0-[x0 y0];
        B=P1-[x0 y0];
        arcL=r*acos(dot(A,B)/(norm(A)*norm(B)));

        L_0_rel(i1,pts(i2))=arcL;
        L_0_rel(pts(i2),i1)=arcL;
        
        
        R_rel(i1,pts(i2))=r;
        R_rel(pts(i2),i1)=r;
                
    end
    
end

