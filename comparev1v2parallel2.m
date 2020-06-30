% Given 2 vectors v1 and v2 which point towards each other, determine how 
% parallel they are. Enfore the condition that v1 and v2 must lie on the
% same half plane in order to be considered

function v1v2p = comparev1v2parallel2(P0,P4,v1,v2)

% First, check that v1 and v2 exist in the same half-plane. If not, make
% v1v2p=1 so that they are not paired together. To do this, construct a
% line from P0 to P4. To do this, identify the endpoints of the vectors v1
% at P0 and v2 at P4
P0t=v1+P0;
P4t=v2+P4;

% determine if each point lies on the same side of the line from P0 to P4
s1=sign((P0t(1)-P0(1))*(P4(2)-P0(2))-(P0t(2)-P0(2))*(P4(1)-P0(1)));
s2=sign((P4t(1)-P0(1))*(P4(2)-P0(2))-(P4t(2)-P0(2))*(P4(1)-P0(1)));
same_side1=s1*s2; % =1 is s1 and s2 have the same sign, i.e., the points are on the same side of the line
% it also =-1 is s1 and s2 are one opposite sides of the line
% it =0 if either s1 or s2 or both are on the line

% apply transformation matrix to move v1 and v2 into quadrant 1 with
% positive vector components
v1_x=v1(1);
v1_y=v1(2);
v2_x=v2(1);
v2_y=v2(2);

T1=[sign(v1_x) 0; 0 sign(v1_y)]; % transformation matrix for v1
T2=[sign(v2_x) 0; 0 sign(v2_y)]; % transformation matrix for v2
v1_temp=v1*T1;
v2_temp=v2*T2;

% STEP 3: Determine how parallel v1_temp and v2_temp are
% normalize each vector
if same_side1>=0 % points are on the same side of the line or at least one point lies on the line
    v1v2p=abs(1-dot(v1_temp./norm(v1_temp),v2_temp./norm(v2_temp)));
else % points are on opposite side of the line
    v1v2p=1;
end


% % CHECKPOINT
% v1v2p
% % construct v1,v2 for plotting
% t=linspace(0,1);
% v1_line_x=t.*(P0(1)+v1(1))+(1-t).*P0(1);
% v1_line_y=t.*(P0(2)+v1(2))+(1-t).*P0(2);
% v2_line_x=t.*(P4(1)+v2(1))+(1-t).*P4(1);
% v2_line_y=t.*(P4(2)+v2(2))+(1-t).*P4(2);
% 
% % construct w1 and w2 for plotting
% w1_1_line_x=t.*(P0(1)+w1_1(1))+(1-t).*P0(1);
% w1_1_line_y=t.*(P0(2)+w1_1(2))+(1-t).*P0(2);
% w1_2_line_x=t.*(P0(1)+w1_2(1))+(1-t).*P0(1);
% w1_2_line_y=t.*(P0(2)+w1_2(2))+(1-t).*P0(2);
% 
% w2_1_line_x=t.*(P4(1)+w2_1(1))+(1-t).*P4(1);
% w2_1_line_y=t.*(P4(2)+w2_1(2))+(1-t).*P4(2);
% w2_2_line_x=t.*(P4(1)+w2_2(1))+(1-t).*P4(1);
% w2_2_line_y=t.*(P4(2)+w2_2(2))+(1-t).*P4(2);
% 
% t=linspace(0,1);
% line_x=t.*(P4(1))+(1-t).*P0(1);
% line_y=t.*(P4(2))+(1-t).*P0(2);

% figure % plot v1,v2, w1s,w2s
% plot(v1_line_x,v1_line_y,'r',v2_line_x,v2_line_y,'r')
% hold on
% plot(P0t(1),P0t(2),'k*',P4t(1),P4t(2),'k*')
% plot(line_x,line_y,'b')
% plot(v1_line_x,v1_line_y,'r',v2_line_x,v2_line_y,'r',w1_1_line_x,w1_1_line_y,'k:',...
%     w1_2_line_x,w1_2_line_y,'k:',w2_1_line_x,w2_1_line_y,'k:',w2_2_line_x,w2_2_line_y,'k:')
% xlim([P0(1),P4(1)])
% ylim([P0(2),P4(2)])
% title('Fixed vectors in red, Variable vectors in black')
% 
% % construct w1i and w2j for plotting
% w1_i_line_x=t.*(P0(1)+w1i(1))+(1-t).*P0(1);
% w1_i_line_y=t.*(P0(2)+w1i(2))+(1-t).*P0(2);
% w2_j_line_x=t.*(P4(1)+w2j(1))+(1-t).*P4(1);
% w2_j_line_y=t.*(P4(2)+w2j(2))+(1-t).*P4(2);
% 
% figure % plot v1,v2, w1i,w2j
% %plot(v1_line_x,v1_line_y,'r')
% plot(v1_line_x,v1_line_y,'r',v2_line_x,v2_line_y,'r',w1_i_line_x,w1_i_line_y,'g:',...
%     w2_j_line_x,w2_j_line_y,'g:')
% %xlim([P0(1),P4(1)])
% %ylim([P0(2),P4(2)])
% title('Fixed vectors in red, "best" vectors in green')
% 
end