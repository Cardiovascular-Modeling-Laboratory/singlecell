function [P1,P2,P3,Fiber_Bending_Energy_std,rtx,rty] = MinEnergyCurveWithoutObstruct3(P0,P4,v1,v2,t,kappa_limit,proceedq,ControlPoints_actin_t,dist_idx,int_idx,vec_ctr)

max_kappa=2*kappa_limit;
%a=7.5;
a1_1=randi([3,10]);%a=8; ;%a=10;
a2_1=randi([3,10]);
%da=5;
max_iters=5;
iter_counter=1;
while max_kappa>kappa_limit && iter_counter<=max_iters
    % compute P1, P3
    if proceedq==1  % if using info from previous time point    
        CPs_old_temp1=ControlPoints_actin_t{dist_idx-1};
        CPs_old_temp2=CPs_old_temp1{int_idx};
        %CPs=CPs_old_temp2{vec_ctr,1};
        CPs=CPs_old_temp2{vec_ctr,1};
        if ~isempty(CPs)
            P1=CPs(2,:);
            P3=CPs(4,:);
        else
            P1 = P0+(norm(P4-P0)/(a1_1*norm(v1))).*v1;
            P3 = P4+(norm(P4-P0)/(a2_1*norm(v2))).*v2;
        end
    else % else, construct P1,P3 from scratch:
        %     P1 = P0+(norm(P4-P0)/(a*norm(v1))).*v1;
        %     P3 = P4+(norm(P4-P0)/(a*norm(v2))).*v2;
        P1 = P0+(norm(P4-P0)/(a1_1*norm(v1))).*v1;
        P3 = P4+(norm(P4-P0)/(a2_1*norm(v2))).*v2;
    end

    % (~Compute P2 and Ebend as follows:)
    Q0 = -P0+2.*P1-2.*P3+P4;
    z1 = 2.*P0-3.*P1-3.*P3+2.*P4;
    z2 = P0+P4;
    a0 = (144/2)*(5/6)*(3/30);
    a1 = (144/2)*(2/3)*(1/30);
    a2 = (144/2)*(1/6)*(5/30);
    
    % Compute P2
    P2x = (a2.*z2(1)-a1.*z1(1))./(2.*(a1+a2));
    P2y = (a2.*z2(2)-a1.*z1(2))./(2.*(a1+a2));
    P2 = [P2x, P2y];

    % Compute Ebend for the standard curve
    Fiber_Bending_Energy_std=a0.*norm(Q0).^2 + a1.*norm(z1+2.*P2).^2 + a2.*norm(z2-2.*P2).^2;
        
    % Construct curve
    rtx = (1-t).^4.*P0(1)+4.*(1-t).^3.*t.*P1(1)+6.*(1-t).^2.*t.^2.*P2(1)+4.*(1-t).*t.^3.*P3(1)+t.^4.*P4(1);
    rty = (1-t).^4.*P0(2)+4.*(1-t).^3.*t.*P1(2)+6.*(1-t).^2.*t.^2.*P2(2)+4.*(1-t).*t.^3.*P3(2)+t.^4.*P4(2);

    CurveLength = arclength(rtx(1,:),rty(1,:));

    
    % Check that curve satisfies max_kappa<= kappa_limit
    Vecs=[rtx' rty'];
    kappa=abs(LineCurvature2D(Vecs));
    max_kappa=max(kappa);
    
    % increment counters and a
%    a=a+da;
    iter_counter=iter_counter+1;
end


end % end function