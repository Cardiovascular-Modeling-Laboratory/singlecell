function combs_new=new_combo_ranking(rank_ij,combo_order)

lin_idx=sub2ind(size(rank_ij),combo_order(:,1),combo_order(:,2));

rank_ij_vals=rank_ij(lin_idx);

[rank_ij_vals_sorted1,sorting_idx1]=sort(rank_ij_vals);
rank_ij_vals_sorted=fliplr(rank_ij_vals_sorted1')';
sorting_idx=fliplr(sorting_idx1')';

combs_new=combo_order(sorting_idx,:);

end % end function