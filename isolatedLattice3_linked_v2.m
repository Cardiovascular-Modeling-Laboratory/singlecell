% This code determines which lattice points need to be removed from
% consideration. It is an extension of the code isolatedLattice.
% First, we determine which lattice points lie within the
% nucleus. Then, we determine which points have no bound integrins, combine
% the indices of these 2 sets into one set and remove the lattice points
% which correspond to those indices.
% User inputs:
%   - An Nx2 lattice, Lat
%   - The radius of the obstruction, r
%   - The x- and y-coordinates of the obstruction, cx and cy
%   - The degree of the lattice points as a row vector, degreeOfPts

function [lat_sat_idx,lat_nuc_idx,lat_other_idx] = isolatedLattice3_linked_v2(Lat,nuc_radius,nuc_cx,nuc_cy,Fsat_lat,V,nuc_rel,R_store,tsec_idx) 
%% Determine which lattice points lie within the nucleus 

dist_sq = (Lat(:,1) -nuc_cx).^2+(Lat(:,2)-nuc_cy).^2;
dist = sqrt(dist_sq);

if nuc_rel==1 % if nucleus should be treated as an obstruction
    lat_nuc_idx = find(dist<=nuc_radius); % lattice numbers of the points which lie within the nucleus
else 
    lat_nuc_idx = [];
end

%% Determine the lattice points which are above the saturation limit
Fmag=sqrt(V(:,1).^2+V(:,2).^2);
r_vals=R_store(max(1,tsec_idx-1),:);
lat_sat_idx=intersect(find(Fmag>=Fsat_lat), find(r_vals>0.2));

%% Determine the remaining lattice points
lat_other_idx_temp=union(lat_nuc_idx,lat_sat_idx);
lat_other_idx=1:size(Lat,1);
lat_other_idx(lat_other_idx_temp)=[];
end