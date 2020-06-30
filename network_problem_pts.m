% determine if there are points with unused vectors that need to be
% adjusted

function pt_remove_vec=network_problem_pts(vecs_used,Lat2,Lat,vec_track,numberOfPts)

pt_remove_vec=[];

% find which elements of vecs_used contain a 1, indiciating that those
% particular vectors in the V2 decomposition have been used and remove them.
% This way, V2_unused will be a cell array containing only those vectors
% which were not used to create a fiber

V2_numbs=cell(1,numberOfPts);
for i=1:numberOfPts
    V2_numbs{i}=1:vec_track(i);
end

for j=1:numel(vecs_used)
    vecs_used_temp=vecs_used{j}; % binary row array indicating which vectors were used for a given point
    used_vec_idx=find(vecs_used_temp==1); % indices of which vectors where used
    vecs_unused_numbs_temp=V2_numbs{j};
    vecs_unused_numbs_temp(used_vec_idx)=[];
    V2_numbs{j}=vecs_unused_numbs_temp;
end


for unused_vec_idx=1:numel(vecs_used)
    vecs_unused_numbs=V2_numbs{unused_vec_idx};
    if ~isempty(vecs_unused_numbs)
        starting_pt=unused_vec_idx; % grab starting saturated lattice point number
        P0 = Lat2(starting_pt,:); % record starting lattice point
        P0_idx2=find(ismember(Lat,P0,'rows')==1);
        pt_remove_vec=unique([pt_remove_vec,P0_idx2]);
    end
end




end %end function