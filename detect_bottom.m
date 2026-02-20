function detect_bottom()

% probably need loop over datagrams 
% might or might not need to check for split ones
% could be do this for single pings only

% what does it need for input



% core of bottom detection
    for aa = 1:Nrx
        if DR(aa) ~= 0
            yBottom(pingidx,aa) = y(aa,DR(aa));
            zBottom(pingidx,aa) = z(aa,DR(aa));
        else
            yBottom(pingidx,aa) = 0;
            zBottom(pingidx,aa) = 0;
        end
    end
    
    
    [~,nadir_idx] = min(abs(yBottom(pingidx,:)-35));
    
    if zBottom(nadir_idx) ~= 0
        nadir_idx = nadir_idx;
    elseif zBottom(nadir_idx + 1) ~= 0
        nadir_idx = nadir_idx + 1;
    else
        nadir_idx = nadir_idx - 1;
    end