function [sim_m, sim_i, m2_sorted] = m2sort(id1, Lat2, V2, m2,theta_ij,outside_segs,dr_dist_squared)

P0=Lat2(id1,:);
W1=V2(id1,:);
sim_m=zeros(1,numel(m2));

for j=1:numel(m2)
    P4=Lat2(m2(j),:);
    W2=V2(m2(j),:);

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
        v_sign=sign(line_eq_m*(P0(1)+W1(:,1)./norm_v)+line_eq_b-(P0(2)+W1(:,2)/norm_v));
        w_sign=sign(line_eq_m*(P4(1)+W2(:,1)./norm_w)+line_eq_b-(P4(2)+W2(:,2)/norm_w));
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
        m=1;
    elseif any(sgn_theta_v.*sgn_theta_w==0)
        sign_idx=find(sgn_theta_v.*sgn_theta_w==0);    
        m=1;
    else
        % this occurs when the vectors for
        % P0 and P4 are seemingly incompatible, i.e, when the curve connecting the
        % 2 points will violate the curvature criterion. to make the code progress,
        % we'll find a vector pair even if the resulting curve won't be used
        sign_idx=find(sgn_theta_v.*sgn_theta_w<0);
        m=NaN;
    end

    % find which w vector is "most similar" to the v vector in the sense that
    % |theta_v-theta_w| is minimal
    sim_m(j)=min(abs(repmat(theta_v,numel(theta_w(sign_idx)),1)-theta_w(sign_idx)))*m;

end

if sum(sum(outside_segs))>0 %if there's a concave portion of the cell
    sim_m_temp=sim_m.*sign(theta_ij(id1,m2));
    sim_i1_temp=find(sim_m_temp<0);
    [~,sim_i1_temp2]=sort(sim_m_temp(sim_i1_temp));
    sim_i1=sim_i1_temp(sim_i1_temp2);
    
    sim_i2_temp=find(sim_m_temp>=0);
    [~,sim_i2_temp2]=sort(sim_m_temp(sim_i2_temp));
    sim_i2=fliplr(sim_i2_temp(sim_i2_temp2));
    
    sim_i=[sim_i1 sim_i2];
    sort_m=[];
%    [~,sim_i]=sort(sim_m.*sign(theta_ij(id1,m2)));
else
    [sort_m,sim_i]=sort(sim_m);
end

if numel(sort_m)~=numel(unique(sort_m))
    sort_m_unique=unique(sort_m,'first');
    for ji=1:numel(sort_m_unique)
        base_idx=find(ismember(sort_m,sort_m_unique(ji))==1);% identify the repeating values
        
        [~,i2]=sort(sqrt(dr_dist_squared(sim_i(base_idx),id1)),'descend'); % sort the repeating values based on the distance between the points, prioritize longer distances over shorter distances

        sim_i(base_idx)=sim_i(base_idx(i2));
    end
end

m2_sorted = m2(sim_i);

end

