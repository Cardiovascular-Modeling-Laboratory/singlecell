% this code determines if the initial fiber constructed hits the nucleus
function [hits_nuc] = avoid_nuc(rtx,rty,nuc_x,nuc_y,nuc_rel)
hits_nuc=0;
if nuc_rel==1 % if nucleus is being treated as an obstruction
    % check if curve hits nucleus
    [x0,y0]=intersections(rtx,rty,nuc_x,nuc_y);
    if ~isempty(x0) || ~isempty(y0)
        hits_nuc=1;
    end      
end

end % end function