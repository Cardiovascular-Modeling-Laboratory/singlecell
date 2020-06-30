function rank_ij = combo_rank(dist_pair,unit_param)

d_ij=dist_pair;
d_ij(d_ij==0)=unit_param;

rank_ij=1-unit_param./d_ij;
end