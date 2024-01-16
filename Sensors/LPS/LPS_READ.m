function [Z_LPS,newLPS,flagLPS] = LPS_READ(LPS_data,mLPS,param)

LPS_dist = LPS_data(mLPS,:);

newLPS = false;
flagLPS = false;
tracking = ~isnan(LPS_dist);
if any(~isnan(LPS_dist))
    newLPS = true;
    flagLPS = true;
end
Z_LPS = {LPS_dist,tracking};

end