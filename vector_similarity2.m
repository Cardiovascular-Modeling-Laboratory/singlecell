% This code determines how similar 2 force vectors are in terms of whether
% they point in the same direction and if they have similar angles. It is
% adapted from wiwjbest3.m
% Given starting/ending points P0,P4 the similarity between vector W2 and
% W1 is determined as
% follows: Draw a line connecting P0 and P4. Determine the angle between W1
% and the vector that points from P0 to P4. For W2, determine the
% angle between that vector and the vector pointing from P4 to P0. The
% angular similarity is determined by taking |theta_v-theta_w|. 
% The direction similarity is determined by using the dot product; if the
% vectors point in the same direction, then the dot product will be 1. if
% they point in different directions, then the dot product will be -1

function [theta_ij_val, delta_ij_val]=vector_similarity2(P0,P4,W1,W2,outline)

% determine if vectors are in the same quadrant, i.e., point in the same
% direction
X1=(W1(1)>=0); 
Y1=(W1(2)>=0);
q1= 1*X1*Y1+2*(1-X1)*Y1+3*(1-X1)*(1-Y1)+4*X1*(1-Y1);

X2=(W2(1)>=0); 
Y2=(W2(2)>=0);
q2= 1*X2*Y2+2*(1-X2)*Y2+3*(1-X2)*(1-Y2)+4*X2*(1-Y2);

if q1==q2
    S_ij_val=1;
else
    S_ij_val=-1;
end

% transform the vectors into quadrant 1 and determine their angular
% similarity
T1=[sign(W1(1)) 0; 0 sign(W1(2))]; % transformation matrix
T2=[sign(W2(1)) 0; 0 sign(W2(2))]; % transformation matrix

W1=W1*T1;
W2=W2*T2;

norm_v=sqrt(W1(:,1).^2+W1(:,2).^2); % norm of v
norm_w=sqrt(W2(:,1).^2+W2(:,2).^2); % norm of w 

%theta_ij_val=real(acos(dot(W1,W2,2)./(norm_v.*norm_w)));
theta_ij_val=S_ij_val*dot(W1,W2,2)./(norm_v.*norm_w);

t=linspace(0,1);
l_x=P0(1).*t+(1-t).*P4(1);
l_y=P0(2).*t+(1-t).*P4(2);

% determine if segment lies on boundary completely. this type of case is
% difficult for inpoly to identify properly so we'll do it manually by
% creating a "sensory region" around the boundary and determining if the 2
% points lie in the same sensory region
pcent=1/100;
P0_segs=zeros(1,size(outline,2)-1);
P4_segs=zeros(1,size(outline,2)-1);
for j=1:size(outline,2)-1
    xmin=min(outline(1,j),outline(1,j+1));
    xmax=max(outline(1,j),outline(1,j+1));
    ymin=min(outline(2,j),outline(2,j+1));
    ymax=max(outline(2,j),outline(2,j+1));
    
    Xmin=xmin*(1-pcent);
    Xmax=xmax*(1+pcent);
    Ymin=ymin*(1-pcent);
    Ymax=ymax*(1+pcent);
    
    P0_segs(j)=Xmin<=P0(1) & P0(1)<=Xmax & Ymin<=P0(2) & P0(2)<=Ymax;
    P4_segs(j)=Xmin<=P4(1) & P4(1)<=Xmax & Ymin<=P4(2) & P4(2)<=Ymax;    
end

if any(P0_segs+P4_segs>1)
    in_determ=1;
else
%     
% 
% % to avoid potential roundoff error that may cause problems if P0,P4 both
% % lie along the same line segment on the boundary, round the line segment
% % to 5 decimal places
% l_x=round(l_x,5);
% l_y=round(l_y,5);
    in_determ=inpoly([l_x' l_y'],round(outline,7)');
    in_determ([1, 100])=[];
end

if any(in_determ==0)
    delta_ij_val=-1;
else
    delta_ij_val=1;
end



% while theta_ij_val>pi
%     theta_ij_val=theta_ij_val-pi;
% end


end % end function