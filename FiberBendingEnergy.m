% this code computes the biological bending energy of a single fiber based
% on the bending energy of the standard bezier curve
function FiberBendingEnergy_curve=FiberBendingEnergy(Fiber_Bending_Energy_std,l_p,CurveLength)

FiberBendingEnergy_curve =(Fiber_Bending_Energy_std)*l_p/(CurveLength)^3 ;