%% this code produces a new way to determine multiple curves after first curve is placed
% Input control points of the known curve and the desired distance between
% the original curve and the new curve to be constructed (t0)

function [P1new, P2new,P3new,Fiber_Bending_Energy_std,rtx,rty] = nudged_curve(P0,P1,P2,P3,P4,t0,t)
%% Determine the 3 step values
%   Get t1,t2,t3 from the values which extremize r' and r''. This means
%   solving sets of polynomials. See notes for where the coefficients of
%   the polynomials came from

a4 = 4.*(P1-P0);
a3 = 12.*(P0-2.*P1+P2);
a2 = -12.*(P0-3.*P1+3.*P2-P3);
a1 = 4.*(P0-4.*P1+6.*P2-4.*P3+P4);

c3 = 12.*(P0-2.*P1+P2);
c2 = -24.*(P0-3.*P1+3.*P2-P3);
c1 = 12.*(P0-4.*P1+6.*P2-4.*P3+P4);

root_set1=roots([a1(1) a2(1) a3(1) a4(1)]);
root_set2=roots([a1(2) a2(2) a3(2) a4(2)]);
root_set3=roots([c1(1) c2(1) c3(1)]);
root_set4=roots([c1(2) c2(2) c3(2)]);
% remove values which are imaginary (so that the final vector can be
% assorted)
root_set1 = root_set1(imag(root_set1)==0);
root_set2 = root_set2(imag(root_set2)==0);
root_set3 = root_set3(imag(root_set3)==0);
root_set4 = root_set4(imag(root_set4)==0);

%   combine all possible solutions into one vector (sorted) and remove
%   values which are outside of (0,1)
root_set = union(union(union(root_set1,root_set2),root_set3),root_set4);
root_set = root_set(root_set>0 & root_set<1);

% round the root values
Ndecimals = 4;
f = 10.^Ndecimals;
root_set = round(f.*root_set)./f;
root_set=unique(root_set);

%   Grab solutions from root_set, depending on how many elements root_set
%   has. Note that root_set is sorted in ascending order
if numel(root_set)==0 % if no elements, predefine ti
    t1=0.25; t2=0.5; t3=0.75;
elseif numel(root_set)==1 % if one element t, other elements are midpoint of [0,t] and [t,1]
    t2=root_set(1);
    t1=t2/2;
    t3=(1+t2)/2;
elseif numel(root_set)==2 % if 2 elements s1,s2, other element is midpoint of [s1,s2]
    t1=root_set(1);
    t3=root_set(2);
    t2=(t1+t3)/2;
elseif numel(root_set)==3 % if 3 elements, done
    t1=root_set(1);
    t2=root_set(2);
    t3=root_set(3);
else % otherwise, there are at least 4 elements. pick smallest, largest, and one close to the median
    t1=root_set(1);
    t3=root_set(end);
    t2 = root_set(ceil(numel(root_set)/2)); 
    % i use this definition because it's equivalent to the median when
    % there are an odd number of elements. when there are an even number of
    % elements, if gives the values that's closest to the median while
    % being smaller than the median
end
t1=0.25; t2=0.5; t3=0.75;
%% Compute the curve
% Construct RHS
%   Compute r1(ti)
r1t1 = ((1-t1)^4).*P0+(4*(1-t1)^3*t1).*P1+(6*(1-t1)^2*t1^2).*P2+(4*(1-t1)*t1^3).*P3+(t1^4).*P4;
r1t2 = ((1-t2)^4).*P0+(4*(1-t2)^3*t2).*P1+(6*(1-t2)^2*t2^2).*P2+(4*(1-t2)*t2^3).*P3+(t2^4).*P4;
r1t3 = ((1-t3)^4).*P0+(4*(1-t3)^3*t3).*P1+(6*(1-t3)^2*t3^2).*P2+(4*(1-t3)*t3^3).*P3+(t3^4).*P4;
%   Compute r1'(ti) and ||r1'(ti)||
r1prime_t1 = (4*(1-t1)^3).*(P1-P0)+(12*(1-t1)^2*t1).*(P2-P1)+(12*(1-t1)*t1^2).*(P3-P2)+(4*t1^3).*(P4-P3);
r1prime_t2 = (4*(1-t2)^3).*(P1-P0)+(12*(1-t2)^2*t2).*(P2-P1)+(12*(1-t2)*t2^2).*(P3-P2)+(4*t2^3).*(P4-P3);
r1prime_t3 = (4*(1-t3)^3).*(P1-P0)+(12*(1-t3)^2*t3).*(P2-P1)+(12*(1-t3)*t3^2).*(P3-P2)+(4*t3^3).*(P4-P3);

r1prime_t1_norm=norm(r1prime_t1);
r1prime_t2_norm=norm(r1prime_t2);
r1prime_t3_norm=norm(r1prime_t3);
%   Compute the normal vectors n(ti)
tangent_t1 = r1prime_t1./(r1prime_t1_norm);
tangent_t2 = r1prime_t2./(r1prime_t2_norm);
tangent_t3 = r1prime_t3./(r1prime_t3_norm);
nt1 = [-tangent_t1(2) tangent_t1(1)];
nt2 = [-tangent_t2(2) tangent_t2(1)];
nt3 = [-tangent_t3(2) tangent_t3(1)];

% here, check that we've got the correct normal vector. if you look at the
% vector P0-P4 and P0-r1t1, then P0-r1t1 is either to the right or left of
% P0-P4. If you look at the normal and tangent vectors, the normal vector
% should be on the same side (right/left) of the tangent vector as P0-r1t1
% is to P0-P4. If not, then change the sign of the normal vector. Check
% that here
% it's possible that the vector P0-P4 is parallel to P0-r1t1. If this is
% the case, then we can nudge in whichever direction we want so long as
% we're consistent

c1test = P4-P0;
c2test1=r1t1-P0;

c2test2=r1t2-P0;
c2test3=r1t3-P0;
s1test1=sign(c1test(1)*(-c2test1(2))+c1test(2)*c2test1(1));
s1test2=sign(c1test(1)*(-c2test2(2))+c1test(2)*c2test2(1));
s1test3=sign(c1test(1)*(-c2test3(2))+c1test(2)*c2test3(1));

s2test1= sign(tangent_t1(1)*(-nt1(2))+tangent_t1(2)*nt1(1));
s2test2= sign(tangent_t2(1)*(-nt2(2))+tangent_t2(2)*nt2(1));
s2test3= sign(tangent_t3(1)*(-nt3(2))+tangent_t3(2)*nt3(1));

paralleltest1=[c1test(2)/c1test(1) c2test1(2)/c2test1(1)];
paralleltest1_val=round(abs(paralleltest1(2)/paralleltest1(1)));
paralleltest2=[c1test(2)/c1test(1) c2test2(2)/c2test2(1)];
paralleltest2_val=round(abs(paralleltest2(2)/paralleltest2(1)));
paralleltest3=[c1test(2)/c1test(1) c2test3(2)/c2test3(1)];
paralleltest3_val=round(abs(paralleltest3(2)/paralleltest3(1)));


if paralleltest1_val==1 && paralleltest2_val==1 && paralleltest3_val==1
    s1=sign(nt1(1));
    s2=sign(nt1(2));
    if sign(nt2(1))~=s1
        nt2(1)=-nt2(1);
    end
    if sign(nt3(1))~=s1
        nt3(1)=-nt3(1);
    end
    if sign(nt2(2))~=s2
        nt2(2)=-nt2(2);
    end
    if sign(nt3(2))~=s2
        nt3(2)=-nt3(2);
    end    
else
    if s1test1~=s2test1
        nt1=-nt1;
    end
    if s1test2~=s2test2
       nt2=-nt2;
    end
    if s1test3~=s2test3
       nt3=-nt3;
    end
end


% if s1test1~=s2test1
%    nt1=-nt1;
% end
% if s1test2~=s2test2
%    nt2=-nt2;
% end
% if s1test3~=s2test3
%    nt3=-nt3;
% end


% if nt1(1)*(-tangent_t1(2))+nt1(2)*tangent_t1(1)>0
%     nt1=-nt1;
% end
% 
% if nt2(1)*(-tangent_t2(2))+nt2(2)*tangent_t2(1)>0
%     nt2=-nt2;
% end
% 
% if nt3(1)*(-tangent_t3(2))+nt3(2)*tangent_t3(1)>0
%     nt3=-nt3;
% end



%   Compute R1(ti)
R1t1= r1t1+t0.*nt1;
R1t2= r1t2+t0.*nt2;
R1t3= r1t3+t0.*nt3;
%   Compute R2(ti)
R2t1 = R1t1-((1-t1)^4).*P0-(t1^4).*P4;
R2t2 = R1t2-((1-t2)^4).*P0-(t2^4).*P4;
R2t3 = R1t3-((1-t3)^4).*P0-(t3^4).*P4;

% Input the components of the coefficient matrix and vectors
f1t1 = 4*(1-t1)^3*t1;
f1t2 = 4*(1-t2)^3*t2;
f1t3 = 4*(1-t3)^3*t3;
f2t1 = 6*(1-t1)^2*t1^2;
f2t2 = 6*(1-t2)^2*t2^2;
f2t3 = 6*(1-t3)^2*t3^2;
f3t1 = 4*(1-t1)*t1^3;
f3t2 = 4*(1-t2)*t2^3;
f3t3 = 4*(1-t3)*t3^3;

R2t1_x= R2t1(1);
R2t2_x= R2t2(1);
R2t3_x= R2t3(1);
R2t1_y= R2t1(2);
R2t2_y= R2t2(2);
R2t3_y= R2t3(2);


% Construct the coefficient matrix and vectors
A = [f1t1 f2t1 f3t1; 
    f1t2 f2t2 f3t2;
    f1t3 f2t3 f3t3];

b1 = [R2t1_x; R2t2_x; R2t3_x];
b2 = [R2t1_y; R2t2_y; R2t3_y];

% record the solutions
X1=linsolve(A,b1);
X2=linsolve(A,b2);

P1new(1)=X1(1);
P2new(1)=X1(2);
P3new(1)=X1(3);
P1new(2)=X2(1);
P2new(2)=X2(2);
P3new(2)=X2(3);

%% Compute bending energy and constructed curve
% (~Compute Ebend as follows:)
Q0 = -P0+2.*P1new-2.*P3new+P4;
z1 = 2.*P0-3.*P1new-3.*P3new+2.*P4;
z2 = P0+P4;
a0 = (144/2)*(5/6)*(3/30);
a1 = (144/2)*(2/3)*(1/30);
a2 = (144/2)*(1/6)*(5/30);
    
% Compute Ebend for the standard curve
Fiber_Bending_Energy_std=a0.*norm(Q0).^2 + a1.*norm(z1+2.*P2new).^2 + a2.*norm(z2-2.*P2new).^2;

% Construct curve
rtx = (1-t).^4.*P0(1)+4.*(1-t).^3.*t.*P1new(1)+6.*(1-t).^2.*t.^2.*P2new(1)+4.*(1-t).*t.^3.*P3new(1)+t.^4.*P4(1);
rty = (1-t).^4.*P0(2)+4.*(1-t).^3.*t.*P1new(2)+6.*(1-t).^2.*t.^2.*P2new(2)+4.*(1-t).*t.^3.*P3new(2)+t.^4.*P4(2);

end