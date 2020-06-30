% this code contains several computations that will need to be computed
% throughout the code. For consistency, they've all been gathered here.
% They include: computation of curve length, computation of the curve's
% bending energy, computation of the maximum curvature in the curve,
% determination if the constructed curve should be kept or not

function [FiberBendingEnergy_curve, keep_curve, curve_prob,max_kappa]=curve_comps(rtx, rty,Fiber_Bending_Energy_std,l_p,P0, P1, P2, P3, P4, t,boundary_pts,kappa_limit)

% compute the curve length
CurveLength = arclength(rtx(1,:),rty(1,:)); 
            
% the fiber bending energy function here is the "normalized"
% bending energy in the sense that it's the bending energy
% (with coefficents l_p*Kb*T) divided by Kb*T
FiberBendingEnergy_curve=FiberBendingEnergy(Fiber_Bending_Energy_std,l_p,CurveLength); % compute fiber bending energy for the biological curve using the bending energy of the standard curve
                    
% compute the biological curvature
Vecs=[rtx' rty'];
%kappa=abs(LineCurvature2D(Vecs)); 
kappa=round(LineCurvature2D(Vecs),5); % record the curvature at each distinct point on the curve and round it 
max_kappa = max(abs(kappa));
                            
% filter computed curve
[keep_curve, curve_prob] = curve_filters(CurveLength,l_p,FiberBendingEnergy_curve,boundary_pts,rtx,rty,kappa_limit,max_kappa,kappa);
end % end function
