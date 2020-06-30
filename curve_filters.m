% this code simply runs through a certain number of filters, i.e., it
% eliminates curves due to some biological criteria. Eliminate curve if:
% - curve if smaller than l_p but has non-zero bending energy
% - curvature of the curve goes above the kappa_limit (this should've
% already been done when constructing curve but check it anyway...)
% - curve goes outside the boundary of the cell
%function [kappa,keep_curve,hits_obstruct ] = curve_filters(CurveLength,l_p,FiberBendingEnergy_curve,boundary_pts,rtx,rty,kappa_limit,kappa,nuc_x,nuc_y,t,fib_sat,Lat,dA_obstruct_radius,fib_sat_limit,P0,P1,P2,P3,P4,bb_area_limit)
function [keep_curve, curve_prob] = curve_filters(CurveLength,l_p,...
    FiberBendingEnergy_curve,boundary_pts,rtx,rty,kappa_limit,max_kappa,kappa)

keep_curve=1; % assume we are keeping the curve
curve_prob=0;

% Filter computed curve through various requirements
% ---> persistence length requirement
if CurveLength<=l_p && FiberBendingEnergy_curve>0 
%    max_kappa=2*kappa_limit;
    keep_curve=0;
    curve_prob=1;
    return
end % end persistence length checkpoint

% ---> curve contained within cell?
% check that the curve (if it exists) doesn't pass outside of the cell. if it does,
% ignore it by setting curvature to be too large. 
% inpoly outputs vector of 1s and 0s (1 if tested point is
% inside polygon, 0 if it's outside polygon). 
test_curve = [rtx' rty'];
in_bdry = inpoly(test_curve,boundary_pts);
% due to possible numerical error, we'll remove the starting and ending
% points of in_bdry from the list. This way, the starting and ending points
% of the curve are not considered
in_bdry([1, 100])=[];



%if any(in_bdry==0)
if sum(abs(in_bdry(1:end-1)-in_bdry(2:end)))==0 % if the curve is completely inside or completely outside cell boundary
%    max_kappa=2*kappa_limit; 
    keep_curve=0;
    curve_prob=2;
%    return
else
    max_kappa=2*kappa_limit;
end % end curve within cell checkpoint

% ---> curvature requirement
ks=sign(kappa);
if max_kappa>kappa_limit || sum(abs(ks(2:end)-ks(1:end-1)))>0
    keep_curve=0;
    curve_prob=3;
    return    
end

end % end function