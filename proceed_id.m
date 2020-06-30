% proceedq: this variable tells the code whether to proceed with the
% curve construction using control points from previous time point (=1) or
% if it should construct curves without using previous time point
% information (=0)

function [proceedq,int_idx] = proceed_id(tsec_idx,Combo_Order_t,starting_pt,ending_pt)

% Update/Set combo order
    if tsec_idx==1 % if this is the first simulation
        proceedq=0; % there's no previous information to use, so set equal to 0
        int_idx=0;
    else % else, pt-pt combos have already been constructed at a previous time point
        combo_order_old=Combo_Order_t{tsec_idx-1};
        
        % check if the pt-pt combo appears in the combo order from the
        % previous time point. Also check if the reverse of the combo order
        % appears. If either is true, then proceedq=1, else proceedq=0
        pt_test=[starting_pt ending_pt; ending_pt starting_pt];
        
        if isempty(combo_order_old) % if the old combo order is empty, then the point shouldn't show up in the previous time point. force this condition:
                combo_order_old=[-1 -1];
        end
            
        [pt_test_result,int_idx]=intersect(combo_order_old,pt_test,'rows');
        
        if ~isempty(pt_test_result) % if pt_test_result is not empty, then this combo shows up in the previous time point
            % use information from previous time to construct the fiber
            proceedq=1;
        else
            proceedq=0;
            int_idx=0;
        end                 
    end
    
end% end function