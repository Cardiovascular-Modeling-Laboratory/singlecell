% This code determine which vector pair to use for a given pair of lattice
% points (P0, P4). W1 is the given vector for P0 and W2 is an array of
% vectors associated with P4.
% Given starting/ending points P0,P4 the appropriate vector w2j from W2
% that should be paired with w1i (equivalent to W1) is determined as
% follows: Draw a line connecting P0 and P4. Determine the angle between W1
% and the vector points from P0 to P4. For each vector in W2, determine the
% angle between that vector and the vector pointing from P4 to P0. The
% vector in W2 that most closely resembles the vector W1 is the vector that
% should be paired (determined by taking |theta_v-theta_w|

function [w1i,w2j, I_col]=wiwjbest3(P0,P4,W1,W2)

% construct the vector pointing from P0 to P4 and the vector pointing from
% P4 to P0
vref=repmat(P4-P0,size(W1,1),1); % P0 to P4: v reference vector
wref=repmat(P0-P4,size(W2,1),1); % P4 to P0: w reference vector

norm_v=sqrt(W1(:,1).^2+W1(:,2).^2); % norm of v
norm_w=sqrt(W2(:,1).^2+W2(:,2).^2); % norm of w vectors

norm_vref=sqrt(vref(:,1).^2+vref(:,2).^2); % norm of v reference vector
norm_wref=sqrt(wref(:,1).^2+wref(:,2).^2); % norm of w reference vector

% determine where P0,P4 are in the plane. this way, vector pairs will only
% be considered if they point towards the same half-plane
if P4(1)-P0(1)~=0
    line_eq_m=(P4(2)-P0(2))/((P4(1)-P0(1)));
    line_eq_b=P4(2)-line_eq_m*P4(1);
    v_sign=sign(line_eq_m*W1(:,1)+line_eq_b-W1(:,2));
    w_sign=sign(line_eq_m*W2(:,1)+line_eq_b-W2(:,2));
else
    v_sign=sign(W1(:,1));
    w_sign=sign(W2(:,1));
end

% compute the angle between v vector and the P0,P4 line; compute the angles
% between each w vector and the P4,P0 line. Note: since i'm using acos,
% it's possible that roundoff error will cause acos to give a complex
% number. so take the real value of the output 
theta_v=real(acos(dot(W1,vref,2)./(norm_v.*norm_vref)).*v_sign);
theta_w=real(acos(dot(W2,wref,2)./(norm_w.*norm_wref)).*w_sign);

% [~,idx_v]=sort(theta_v);
% [W1(idx_v,:),theta_v(idx_v)]
% 
% [~,idx_w]=sort(theta_w);
% [W2(idx_w,:),theta_w(idx_w)]

sgn_theta_v=repmat(sign(theta_v),numel(theta_w),1);
sgn_theta_w=sign(theta_w);

% isolate w vectors that are in the same half-plane as v
if any(sgn_theta_v.*sgn_theta_w>0)
    %sign_idx=find(sign(theta_w)==sign(theta_v)); 
    sign_idx=find(sgn_theta_v.*sgn_theta_w>0);
elseif any(sgn_theta_v.*sgn_theta_w==0)
    sign_idx=find(sgn_theta_v.*sgn_theta_w==0);    
else
    % this occurs when the vectors for
    % P0 and P4 are seemingly incompatible, i.e, when the curve connecting the
    % 2 points will violate the curvature criterion. to make the code progress,
    % we'll find a vector pair even if the resulting curve won't be used
    sign_idx=find(sgn_theta_v.*sgn_theta_w<0);
end

% find which w vector is "most similar" to the v vector in the sense that
% |theta_v-theta_w| is minimal
[~,min_idx]=min(abs(repmat(theta_v,numel(theta_w(sign_idx)),1)-theta_w(sign_idx)));

% grab the desired vector from the W2 array
w2j=W2(sign_idx(min_idx(1)),:);
I_col = find(ismember(W2,w2j,'rows')==1); % record which row of W2 contained w2j (don't pay attention to the labeling as "I_col")
w1i=W1;

end % end function