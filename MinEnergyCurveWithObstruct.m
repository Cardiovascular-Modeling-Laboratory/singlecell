function [P2_1,P2_2,Fiber_Bending_Energy_std_1,Fiber_Bending_Energy_std_2,rtx_1,rty_1,rtx_2,rty_2] = MinEnergyCurveWithObstruct(P0,P1,P3,P4,t,cx,cy,radius)
ctr=[cx, cy];

Q0 = -P0+2.*P1-2.*P3+P4;
z1 = 2.*P0-3.*P1-3.*P3+2.*P4;
z2 = P0+P4;
a0 = (144/2)*(5/6)*(3/30);
a1 = (144/2)*(2/3)*(1/30);
a2 = (144/2)*(1/6)*(5/30);

%% Compute points of intersection with obstruction
Pc1=zeros(1,2);
Pc2=zeros(1,2);

% check that P0 & P4 don't lie on a horizonal line
if P0(2)-P4(2)==0
    Pc1(1) = ctr(1);
    Pc2(1) = ctr(1);
    % for x, there are 2 possible y's
    Pc1(2) = ctr(2)+radius;
    Pc2(2)= ctr(2)-radius;
else
    m = (P4(2)-P0(2))/(P4(1)-P0(1));
    b2 = ctr(2)+(1/m)*ctr(1);
    A = 1+1/m^2;
    B = -((2/m)*(b2-ctr(2))+2*ctr(1));
    C = (b2-ctr(2))^2+ctr(1)^2-radius^2;
    pt2_x = roots([A B C]); % find roots of quadratic 
    pt2_x = pt2_x(imag(pt2_x)==0); % isolate only real roots, there will be 2
    Pc1(1) = pt2_x(1); % root #1
    Pc2(1) = pt2_x(2); % root #2

    Pc1(2) = -(1/m)*Pc1(1)+b2; % y value associated with root #1
    Pc2(2) = -(1/m)*Pc2(1)+b2; % y value associated with root #2
end % end if statement for computing intersection points

%% Compute T value which minimizes bending energy for curves 1 and 2
% 		Compute time point solutions (can be up to 8 time points) for
% 		each Pc
        q = a1.*z1-a2.*z2;
% polynomial coefficients for Pc1:
    p11=(-a1-a2)*dot(P0,P0)+6*(a1+a2)*dot(P0,P1)+2*(a1+a2)*dot(P0,P3)-8*(a1+a2)*dot(P1,P1)...
           -2*(a1+a2)*dot(P1,P4)+8*(a1+a2)*dot(P3,P3)-6*(a1+a2)*dot(P3,P4)+(a1+a2)*dot(P4,P4)...
           +3*dot(q,P0)-6*dot(q,P1)+6*dot(q,P3)-3*dot(q,P4);

       p21=8*(a1+a2)*dot(P0,P0)-42*(a1+a2)*dot(P0,P1)-10*(a1+a2)*dot(P0,P3)+48*(a1+a2)*dot(P1,P1)...
           +6*(a1+a2)*dot(P1,P4)-16*(a1+a2)*dot(P3,P3)+6*(a1+a2)*dot(P3,P4)-18*dot(q,P0)...
           +30*dot(q,P1)-18*dot(q,P3)+6*dot(q,P4);

       p31=-28*(a1+a2)*dot(P0,P0)+126*(a1+a2)*dot(P0,P1)+20*(a1+a2)*dot(P0,P3)...
           -120*(a1+a2)*dot(P1,P1)-6*(a1+a2)*dot(P1,P4)+8*(a1+a2)*dot(P3,P3)...
           +45*dot(q,P0)-60*dot(q,P1)+18*dot(q,P3)-3*dot(q,P4);

       p41=56*(a1+a2)*dot(P0,P0)-210*(a1+a2)*dot(P0,P1)-20*(a1+a2)*dot(P0,P3)-2*(a1+a2)*dot(P0,Pc1)...
           +160*(a1+a2)*dot(P1,P1)+2*(a1+a2)*dot(P1,P4)+8*(a1+a2)*dot(P1,Pc1)+8*(a1+a2)*dot(P3,Pc1)...
           -2*(a1+a2)*dot(P4,Pc1)-60*dot(q,P0)+60*dot(q,P1)-6*dot(q,P3)+6*dot(q,Pc1);

       p51= -70*(a1+a2)*dot(P0,P0)+210*(a1+a2)*dot(P0,P1)+10*(a1+a2)*dot(P0,P3)...
           +10*(a1+a2)*dot(P0,Pc1)-120*(a1+a2)*dot(P1,P1)-30*(a1+a2)*dot(P1,Pc1)...
           -10*(a1+a2)*dot(P3,Pc1)+45*dot(q,P0)-30*dot(q,P1)-15*dot(q,Pc1);

       p61=56*(a1+a2)*dot(P0,P0)-126*(a1+a2)*dot(P0,P1)-2*(a1+a2)*dot(P0,P3)-20*(a1+a2)*dot(P0,Pc1)...
           +48*(a1+a2)*dot(P1,P1)+42*(a1+a2)*dot(P1,Pc1)+2*(a1+a2)*dot(P3,Pc1)-18*dot(q,P0)...
           +6*dot(q,P1)+12*dot(q,Pc1);

       p71=-28*(a1+a2)*dot(P0,P0)+42*(a1+a2)*dot(P0,P1)+20*(a1+a2)*dot(P0,Pc1)-8*(a1+a2)*dot(P1,P1)...
           -26*(a1+a2)*dot(P1,Pc1)+3*dot(q,P0)-3*dot(q,Pc1);

       p81=8*(a1+a2)*dot(P0,P0)-6*(a1+a2)*dot(P0,P1)-10*(a1+a2)*dot(P0,Pc1)...
       +6*(a1+a2)*dot(P1,Pc1)+2*(a1+a2)*dot(Pc1,Pc1);

        p91=-(a1+a2)*dot(P0,P0)+2*(a1+a2)*dot(P0,Pc1)-(a1+a2)*dot(Pc1,Pc1);

% polynomial coefficients for Pc2:
       p12=(-a1-a2)*dot(P0,P0)+6*(a1+a2)*dot(P0,P1)+2*(a1+a2)*dot(P0,P3)-8*(a1+a2)*dot(P1,P1)...
           -2*(a1+a2)*dot(P1,P4)+8*(a1+a2)*dot(P3,P3)-6*(a1+a2)*dot(P3,P4)+(a1+a2)*dot(P4,P4)...
           +3*dot(q,P0)-6*dot(q,P1)+6*dot(q,P3)-3*dot(q,P4);

       p22=8*(a1+a2)*dot(P0,P0)-42*(a1+a2)*dot(P0,P1)-10*(a1+a2)*dot(P0,P3)+48*(a1+a2)*dot(P1,P1)...
           +6*(a1+a2)*dot(P1,P4)-16*(a1+a2)*dot(P3,P3)+6*(a1+a2)*dot(P3,P4)-18*dot(q,P0)...
           +30*dot(q,P1)-18*dot(q,P3)+6*dot(q,P4);

       p32=-28*(a1+a2)*dot(P0,P0)+126*(a1+a2)*dot(P0,P1)+20*(a1+a2)*dot(P0,P3)...
           -120*(a1+a2)*dot(P1,P1)-6*(a1+a2)*dot(P1,P4)+8*(a1+a2)*dot(P3,P3)...
           +45*dot(q,P0)-60*dot(q,P1)+18*dot(q,P3)-3*dot(q,P4);

       p42=56*(a1+a2)*dot(P0,P0)-210*(a1+a2)*dot(P0,P1)-20*(a1+a2)*dot(P0,P3)-2*(a1+a2)*dot(P0,Pc2)...
           +160*(a1+a2)*dot(P1,P1)+2*(a1+a2)*dot(P1,P4)+8*(a1+a2)*dot(P1,Pc2)+8*(a1+a2)*dot(P3,Pc2)...
           -2*(a1+a2)*dot(P4,Pc2)-60*dot(q,P0)+60*dot(q,P1)-6*dot(q,P3)+6*dot(q,Pc2);

       p52= -70*(a1+a2)*dot(P0,P0)+210*(a1+a2)*dot(P0,P1)+10*(a1+a2)*dot(P0,P3)...
           +10*(a1+a2)*dot(P0,Pc2)-120*(a1+a2)*dot(P1,P1)-30*(a1+a2)*dot(P1,Pc2)...
           -10*(a1+a2)*dot(P3,Pc2)+45*dot(q,P0)-30*dot(q,P1)-15*dot(q,Pc2);

       p62=56*(a1+a2)*dot(P0,P0)-126*(a1+a2)*dot(P0,P1)-2*(a1+a2)*dot(P0,P3)-20*(a1+a2)*dot(P0,Pc2)...
           +48*(a1+a2)*dot(P1,P1)+42*(a1+a2)*dot(P1,Pc2)+2*(a1+a2)*dot(P3,Pc2)-18*dot(q,P0)...
           +6*dot(q,P1)+12*dot(q,Pc2);

       p72=-28*(a1+a2)*dot(P0,P0)+42*(a1+a2)*dot(P0,P1)+20*(a1+a2)*dot(P0,Pc2)-8*(a1+a2)*dot(P1,P1)...
           -26*(a1+a2)*dot(P1,Pc2)+3*dot(q,P0)-3*dot(q,Pc2);

       p82=8*(a1+a2)*dot(P0,P0)-6*(a1+a2)*dot(P0,P1)-10*(a1+a2)*dot(P0,Pc2)...
       +6*(a1+a2)*dot(P1,Pc2)+2*(a1+a2)*dot(Pc2,Pc2);

        p92=-(a1+a2)*dot(P0,P0)+2*(a1+a2)*dot(P0,Pc2)-(a1+a2)*dot(Pc2,Pc2);

% find roots for each polynomial (will be displayed as a column vector if there are multiple solutions):    
T_Pc1=roots([p11, p21, p31, p41, p51, p61, p71, p81, p91]);
T_Pc2=roots([p12, p22, p32, p42, p52, p62, p72, p82, p92]);

% 		Isolate time point solutions that are Real and in [0,1]
T_Pc1 = T_Pc1(imag(T_Pc1)==0);
T_Pc1 = T_Pc1(T_Pc1>=0 & T_Pc1<=1);
T_Pc2 = T_Pc2(imag(T_Pc2)==0);
T_Pc2 = T_Pc2(T_Pc2>=0 & T_Pc2<=1);

if isempty(T_Pc1)
    disp('error in computing min energy curve with obstruct')
    T_Pc1 = roots([p11, p21, p31, p41, p51, p61, p71, p81, p91])
    T_Pc1 = T_Pc1(imag(T_Pc1)==0)
end

%% compute P2_1
% construct P2_1 for each associated step T_Pc1 (there can be up to 8)
P2_1=zeros(length(T_Pc1),2);
for i=1:length(T_Pc1)
    T=T_Pc1(i);
    P2_1(i,:) = (Pc1-(1-T).^4.*P0-4.*(1-T).^3.*T.*P1-4.*(1-T).*T.^3.*P3-T^4.*P4)./(6.*(1-T).^2.*T.^2);
end

% 		Compute Ebend_1 for P2_1
Ebend_1 = zeros(1,length(T_Pc1));    
for i=1:length(Ebend_1)
    Ebend_1(i)=a0.*norm(Q0).^2 + a1.*norm(z1+2.*P2_1(i,:)).^2 + a2.*norm(z2-2.*P2_1(i,:)).^2;
end

% Since we're using approximate values, round all the bending energies
% to 8 decimal places
Ndecimals = 8;
f = 10.^Ndecimals;
Ebend_1 = round(f.*Ebend_1)./f;

% 		Find P2_1 which gives min(Ebend) 
idx=find(Ebend_1==min(Ebend_1));
P2_1=P2_1(idx(1),:);
Ebend_1=min(Ebend_1); %



%% compute P2_2
% construct P2_2 for each associated step T_Pc2 (there can be up to 8)
P2_2=zeros(length(T_Pc2),2);
for i=1:length(T_Pc2)
    T=T_Pc2(i);
    P2_2(i,:) = (Pc2-(1-T).^4.*P0-4.*(1-T).^3.*T.*P1-4.*(1-T).*T.^3.*P3-T^4.*P4)./(6.*(1-T).^2.*T.^2);
end

% 		Compute Ebend_2 for P2_2
Ebend_2 = zeros(1,length(T_Pc2));    
for i=1:length(Ebend_2)
    Ebend_2(i)=a0.*norm(Q0).^2 + a1.*norm(z1+2.*P2_2(i,:)).^2 + a2.*norm(z2-2.*P2_2(i,:)).^2;
end

% Since we're using approximate values, round all the bending energies
% to 8 decimal places
Ndecimals = 8;
f = 10.^Ndecimals;
Ebend_2 = round(f.*Ebend_2)./f;

% 		Find P2_2 which gives min(Ebend) 
idx=find(Ebend_2==min(Ebend_2));
P2_2=P2_2(idx(1),:);
Ebend_2=min(Ebend_2); %


%% construct curve 1 

% Compute Ebend for the standard curve (curve 1)
Fiber_Bending_Energy_std_1=Ebend_1;

% Construct curve 1
rtx_1 = (1-t).^4.*P0(1)+4.*(1-t).^3.*t.*P1(1)+6.*(1-t).^2.*t.^2.*P2_1(1)+4.*(1-t).*t.^3.*P3(1)+t.^4.*P4(1);
rty_1 = (1-t).^4.*P0(2)+4.*(1-t).^3.*t.*P1(2)+6.*(1-t).^2.*t.^2.*P2_1(2)+4.*(1-t).*t.^3.*P3(2)+t.^4.*P4(2);

%% construct curve 2

% Compute Ebend for the standard curve (curve 2)
Fiber_Bending_Energy_std_2=Ebend_2;

% Construct curve 2
rtx_2 = (1-t).^4.*P0(1)+4.*(1-t).^3.*t.*P1(1)+6.*(1-t).^2.*t.^2.*P2_2(1)+4.*(1-t).*t.^3.*P3(1)+t.^4.*P4(1);
rty_2 = (1-t).^4.*P0(2)+4.*(1-t).^3.*t.*P1(2)+6.*(1-t).^2.*t.^2.*P2_2(2)+4.*(1-t).*t.^3.*P3(2)+t.^4.*P4(2);



end