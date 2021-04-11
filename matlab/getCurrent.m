function [current,ori,pnti,fwd] = getCurrent(pow,pnt,M,gain, maxdist, threshold)
%   Input:
% pow is locations by 1 vector with dipole strength
% pnt, M and gain are the same as in rimda.m
% recomended:
% maxdist dist is maximum distance for local maxima search (30 mm);
% threshold is the ratio of the maximum dipole, below which other dipoles are ignored (0.3);
%   Output:
% pnti, current and fwd are the same as in rimda.m
% ori is moment orientation
%% find local maxima
% selsct local maxima in space
Pow2 = sqrt(pow(1:length(pnt),:).^2+pow(length(pnt)+1:length(pnt)*2,:).^2);
maxima = false(size(Pow2));
for pnti = 1:length(maxima)
    pos = pnt(pnti,:);
    distnc = sqrt(sum((pnt-repmat(pos,length(pnt),1)).^2,2));
    neighb = distnc<maxdist;
    for sampi = 1:size(M,2)
        powNeighb = max(Pow2(neighb,sampi));
        if sampi>1
            powNeighb = [powNeighb,max(Pow2(neighb,sampi-1))];
        end
        if sampi<size(M,2)
            powNeighb = [powNeighb,max(Pow2(neighb,sampi+1))];
        end
        maxPow(sampi) = max(powNeighb);
    end
    maxTimes = Pow2(pnti,:)==maxPow;
    maxima(pnti,maxTimes) = true;
end
maxima1 = Pow2.*maxima;
maxima2 = maxima1>(max(max(maxima1))*threshold);
pnti = find(sum(maxima2,2)>0);

% find timepoint
times = zeros(size(pnti));
for maxj = 1:length(pnti)
    maxm = find(maxima2(pnti(maxj),:));
    [~,maxl] = max(maxima1(pnti(maxj),maxm));
    times(maxj) = maxm(maxl);
end
Gain = [];
for maxi = 1:length(pnti)
    maxj = pnti(maxi);
    maxk = pnti(maxi)+length(Pow2);
    normFac = sqrt(1/((pow(maxj,times(maxi))).^2+(pow(maxk,times(maxi))).^2));
    ori(maxi,1:2) = [pow(maxj,times(maxi))*normFac pow(maxk,times(maxi))*normFac];
    Gain(1:248,maxi) = gain(:,maxj)*ori(maxi,1)+gain(:,maxk)*ori(maxi,2);
end
current = Gain\M;
fwd = Gain*current;
R = (corr(fwd(:),M(:))).^2;


