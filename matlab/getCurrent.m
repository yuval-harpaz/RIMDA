function [current,ori,pnti,fwd]=getCurrent2(pow,pnt,M,gain, maxdist, threshold, figs)

% pow is Pow1, 2 dipoles per location
% recomended:
% maxdist=30;
% threshold=0.5;
% figs can be true or false
% current=[];
sampBefore=1;
%% find local maxima
if ~exist('figs','var')
    figs=true;
end
% % selsct local maxima in space
% Pow2=sqrt(pow(1:length(pnt),:).^2+pow(length(pnt)+1:length(pnt)*2,:).^2);
% maxima=false(size(Pow2));
% for pnti=1:length(maxima)
%     pos=pnt(pnti,:);
%     distnc=sqrt(sum((pnt-repmat(pos,length(pnt),1)).^2,2));
%     neighb=distnc<maxdist;
%     maxTimes=Pow2(pnti,:)==max(Pow2(neighb,:));
%     maxima(pnti,maxTimes)=true;
% end


% selsct local maxima in space
Pow2=sqrt(pow(1:length(pnt),:).^2+pow(length(pnt)+1:length(pnt)*2,:).^2);
maxima=false(size(Pow2));
for pnti=1:length(maxima)
    pos=pnt(pnti,:);
    distnc=sqrt(sum((pnt-repmat(pos,length(pnt),1)).^2,2));
    neighb=distnc<maxdist;
    for sampi=1:size(M,2)
        powNeighb=max(Pow2(neighb,sampi));
        if sampi>1
            powNeighb=[powNeighb,max(Pow2(neighb,sampi-1))];
        end
        if sampi<size(M,2)
            powNeighb=[powNeighb,max(Pow2(neighb,sampi+1))];
        end
        maxPow(sampi)=max(powNeighb);
    end
    maxTimes=Pow2(pnti,:)==maxPow;
    maxima(pnti,maxTimes)=true;
end

maxima1=Pow2.*maxima;
%maxima2=maxima1>repmat(max(maxima1),size(maxima1,1),1)*threshold;
maxima2=maxima1>(max(max(maxima1))*threshold);
pnti=find(sum(maxima2,2)>0);

% find timepoint
times=zeros(size(pnti));
for maxj=1:length(pnti)
    maxm=find(maxima2(pnti(maxj),:));
    [~,maxl]=max(maxima1(pnti(maxj),maxm));
    times(maxj)=maxm(maxl);
end

% [~,maxOrder]=sort(Pow2(maxima),'descend');
% maxima=maxima(maxOrder);
% 
% pp=Pow2(maxima);
% pp=pp/max(pp);
% if figs
%     disp([' number of local maxima found ', num2str(size(maxima,1))]);
%     disp(pp)
% end
% pnti=maxima(pp>=threshold);

% right is 777, left is 757
Gain=[];
for maxi=1:length(pnti)
    maxj=pnti(maxi);
    maxk=pnti(maxi)+length(Pow2);
    normFac=sqrt(1/((pow(maxj,times(maxi))).^2+(pow(maxk,times(maxi))).^2));
    ori(maxi,1:2)=[pow(maxj,times(maxi))*normFac pow(maxk,times(maxi))*normFac];
    Gain(1:248,maxi)=gain(:,maxj)*ori(maxi,1)+gain(:,maxk)*ori(maxi,2);
end
%source=Gain\M;
% R=[];
% for maxi=1:length(maxima)
%     source=Gain(:,1:maxi)\M;
%     R(maxi)=(corr(Gain(:,1:maxi)*source,M)).^2;
% end
% R
% figure;
% scatter3pnt(pnt,25,Pow2)
current=Gain\M;
fwd=Gain*current;
R=(corr(fwd(:),M(:))).^2;


% if figs
%     src=zeros(size(Pow2,1),1);
%     src(pnti)=mean(abs(current),2);
%     figure;
%     scatter3pnt(pnt,25,src)
%     figure;
%     if size(fwd,2)==1
%         topoplot248(fwd);
%     else
%         topoplot248(mean(abs(fwd),2));
%     end
%     title([num2str(size(pnti,1)),' dipoles explain ',num2str(round(R*100)),'%'])
% end
