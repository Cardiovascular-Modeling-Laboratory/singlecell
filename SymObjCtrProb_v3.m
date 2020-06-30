% This function determines the center of one of the 2D symmetric geometries
% used cellgeometry code. 
% Input the minimum number of lattice points, the area of the cell,
% the "choice" number to determine cell geometry, the length of the
% major/minor axes of the confidence ellipse region, the radius of the 
% obstruction, and the angle (in degrees) the the major axis of the
% confidence regions makes with the x-axis, and the geometry info (lattice
% points and shape choice).
% The code determines the exact center of the geometry. Then, for
% variability, the center is considered to be the average of a 2d gaussian
% distribution. A random point is then picked from the probability 
% distribution. If the point falls outside of the geometry, pick a 
% different point. If the point lies within the geometry, call that point 
% the center of the object. 
% As of 3/24/2017, the code is capable of determining the center for the 
% following geometries: 
% Square/Rectangle, Circle/Oval, Triangle
% Note: this code requires the cell to be in the micrometer scale so
% convert it

function [nuc_cx, nuc_cy] = SymObjCtrProb_v3(A,choice,nuc_radius,mat_r)
A=A*10^(12); % convert to micometer scale
nuc_radius=nuc_radius*10^6;
mat_r=mat_r.*10^6;

% base parameter set
r1=2; % length of major axis of confidence ellipse (in micrometer)
r2=2; % length of minor axis of confidence ellipse (in micrometer)
angle=0;  % angle of major axis of confidence ellipse 

% ----------------- Construct sigma
% Construct sigma so that 68% of randomly selected points lie within a
% specified ellipse (or circle if r1=r2)
X2=2.28; %68% confidence, 2 degrees of freedom: 1-alpha=0.68 --> alpha=0.32 --> X2=2.28, obtained from http://www.mathcelebrity.com/chicritical.php?a=0.99&d=14&pl=Calculate
%X2=5.99; %95% confidence, 2 degrees of freedom: 1-alpha=0.95 --> alpha=0.05 --> X2=5.99
%X2=3.22; %80% confidence, 2 degrees of freedom: 1-alpha=0.8 --> alpha=0.2 --> X2=3.22

if angle==0
    sigma = [r1^2/X2 0; 0 r2^2/X2] ;
elseif angle==90 || angle==-90
    sigma = [r2^2/X2 0; 0 r1^2/X2];
else
    l1=(r1*cosd(angle))^2/X2;
    l2=(r2*cosd(angle))^2/X2;
    sigma = [(l1+l2*(tand(angle)^2))/(secd(angle)^2) (l1-l2)*tand(angle)/(secd(angle)^2) ...
            (l1-l2)*tand(angle)/(secd(angle)^2) (l2+l1*(tand(angle)^2))/(secd(angle)^2)];
end

k=1; % index used to stop point testing from being an infinite while loop

xyMin=min(mat_r); % vector contain xmin and ymin
xyMax=max(mat_r); % vector contain xmax and ymax

xmin=xyMin(1); % isolate xmin
xmax=xyMax(1); % isolate xmax
ymin=xyMin(2); % isolate ymin
ymax=xyMax(2); % isolate ymax
    
if choice==1 || choice ==2 || choice==4% if geometry is a square, rectangle, triangle
    % Determine exact center of geometry
    l=xmax-xmin;
    w=ymax-ymin;

    ExactCtr = [xmin+l/2, ymin+w/2]; % the exact center of the geometry
    
    % Pick point from 2d gaussian distribution averaged at the exact center
    mu = ExactCtr;
    
    % Pick a point at random from the distribution
    ctr = mvnrnd(mu,sigma,1);
    
    % Determine if the obstruction lies within the geometry
    if choice==1 || choice==2 % if in a square or rectangle...
        % determine if the left/right/top/bottom of nuclues are within the
        % left/right/top/bottom of geometry
            nucleus_xmin=ctr(1)-nuc_radius;
            nucleus_xmax=ctr(1)+nuc_radius;
            nucleus_ymin=ctr(2)-nuc_radius;
            nucleus_ymax=ctr(2)+nuc_radius;
        while nucleus_xmin<xmin || nucleus_xmax>xmax || nucleus_ymin<ymin || nucleus_ymax>ymax
            if k>15 % attempted 15 times to place nucleus, kept getting errors so just place it in the center
                ctr=ExactCtr;
                %disp('Error: Made 15 attempts to place the nucleus but no placement lies within the square/rectangle');
                break
            end
            k=k+1;
            ctr = mvnrnd(mu,sigma,1);
        end
    elseif choice==4 % if in a triangle
        % rotate the left/right side of the triangle and the center of the
        % nucleus by the appropriate angles
        ctrx=ctr(1);
        ctry=ctr(2);
        idx=find(mat_r(:,2)==max(mat_r(:,2)));
        x1=mat_r(idx,1);
        
        % rotation of left side of triangle clockwise by angle created
        th1=atand((ymax-ymin)/(x1-xmin));
        R1=[cos(th1) sin(th1); -sin(th1) cos(th1)];
        ctr_tilde1=R1*[ctrx ctry]'; % new rotated center
        
        % rotation of right side of triangle counterclockwise by angle
        % created
        th2=atand((ymax-ymin)/(x1-xmax));
        R2=[cos(th2) -sin(th2); sin(th2) cos(th2)];
        ctr_tilde2=R2*[ctrx ctry]'; % new rotated center
        
        % determine if nucleus is above the bottom line, to the right of
        % the left line (using new rotated coordinates), to the left of the
        % right line (using new rotated coordinates)
        while ctry-ymin-nuc_radius<0 || ctr_tilde1(2)-ymin+nuc_radius>0 || ctr_tilde2(2)-ymin+nuc_radius>0
            if k>15
                ctr=ExactCtr;
%                disp('Error: Made 15 attempts to place the nucleus but no placement lies within the triangle');
                break
            end
            k=k+1;
            ctr = mvnrnd(mu,sigma,1);
        end
    end
elseif choice==3 || choice==7
    % Determine exact center of geometry
    ExactCtr = [0 0]; % exact center of the geometry 
    
    % Pick point from 2d gaussian distribution averaged at the exact center
    mu = ExactCtr;

    % Pick a point at random from the distribution
    ctr = mvnrnd(mu,sigma,1);

	nuc_cx=ctr(1);
	nuc_cy=ctr(2);
	% intersection point 1 of the line thru center of main geometry and
	% center of nucleus
    x1= (nuc_cx^3+nuc_cx*nuc_cy^2-sqrt( nuc_cx^2*(nuc_cx^2+nuc_cy^2)*nuc_radius^2 ))/(nuc_cx^2+nuc_cy^2);
    y1= (nuc_cx^3*nuc_cy+nuc_cx*nuc_cy^3-nuc_cy*sqrt( nuc_cx^2*(nuc_cx^2+nuc_cy^2)*nuc_radius^2 ))/(nuc_cx^2+nuc_cy^2);

    % intersection point 2 of the line thru center of main geometry and
    % center of nucleus
    x2= (nuc_cx^3+nuc_cx*nuc_cy^2+sqrt(nuc_cx^2*(nuc_cx^2+nuc_cy^2)*nuc_radius^2))/(nuc_cx^2+nuc_cy^2);
    y2= (nuc_cy*(nuc_cx^3+nuc_cx*nuc_cy^2+sqrt(nuc_cx^2*(nuc_cx^2+nuc_cy^2)*nuc_radius^2)))/(nuc_cx*(nuc_cx^2+nuc_cy^2));
        
    % distance from center of main geometry to intersection points.
    % Pick the maximum of these 2 distances
    d1=sqrt(x1^2+y1^2);
    d2=sqrt(x2^2+y2^2);
    dmax=max(d1,d2);

    
    % Determine if the obstruction lies within the geometry
    if choice==3 % if in a circle
        r=1; % radius of the main cirular geometry. This should be automated in later versions

        while dmax>r%temp condition
            if k>15
                ctr=ExactCtr;
%                disp('Error: Made 15 attempts to place the nucleus but no placement lies within the circle');
                break
            end
            k=k+1;
            ctr = mvnrnd(mu,sigma,1);
        end    
    elseif choice==7 % if in an ellipse
        Rx=2; %major axis of ellipse geometry. this should be automated in later versions
        Ry=1; %minor axis of ellipse geometry. this should be automated in later versions
        
        % Determine where the line connecting the center of main geometry
        % and the center of the nucleus intersects the ellipse
        x3 = -(nuc_cx*Rx*Ry)/sqrt( nuc_cy^2*Rx^2+nuc_cx^2*Ry^2 );
        y3 = -(nuc_cy*Rx*Ry)/sqrt( nuc_cy^2*Rx^2+nuc_cx^2*Ry^2 );
        
        x4 = -x3;
        y4 = -y3;
        
        % Determine which of these two points lies on the line segment that
        % goes through the center of the nucleus. To do this, the sign of
        % the points that intersect with the ellipse should have the same
        % sign as the points that intersect with the center of the nucleus
        
        A3 = [x3 y3; x4 y4];
        b=[nuc_cx nuc_cy];
        id = strfind(reshape(sign(A)',1,[],sign(b)));
        interPt = A3(id,:);
        
        % Determine the distance between the center of the main geometry
        % and this intersection point.
        d_tilde = sqrt(interPt(1)^2+interPt(2)^2);
        
        while dmax>d_tilde %condition
            if k>15
                ctr=ExactCtr;
%                disp('Error: Made 15 attempts to place the nucleus but no placement lies within the ellipse');
                break
            end
            k=k+1;
            ctr = mvnrnd(mu,sigma,1);
        end    
    end
end

nuc_cx=ctr(1);
nuc_cy=ctr(2);

% convert center back to meter scale
nuc_cx=nuc_cx*10^(-6);
nuc_cy=nuc_cy*10^(-6);

% % plot the nucleus, and the resulting 2d gaussian distribution
% figure
% %subplot(1,2,1)
%     t=linspace(0,2*pi);
%     xt = radius*cos(t)+ctr_x;
%     yt = radius*sin(t)+ctr_y;
%     %plot(mat_r(:,1),mat_r(:,2),'r*',ExactCtr(1), ExactCtr(2),'b*') % plot lattice and the exact center of the lattice
%     plot(ExactCtr(1), ExactCtr(2),'b*') % plot the exact center of the geometry
%     hold on
%     %plot(xt, yt, 'g') % plot nucleus
%     patch(xt,yt,'blue','facealpha',0.3);
%     axis equal
% set(gca,'xtick',[])
% set(gca,'xticklabel',[])
% set(gca,'ytick',[])
% set(gca,'yticklabel',[])
% 
% tmp = size(outline);
% hold on;
% if tmp(1) ==1 %if circle
%     rectangle('Position',[outline(1)-outline(3),outline(2)-outline(3),2*outline(3),2*outline(3)],'Curvature',[1,1]);%,...
%     %          'FaceColor','r')
% else
%     line(outline(1,:),outline(2,:),'Color',[0 0 0],'LineWidth',2);
% end    


end %function end



