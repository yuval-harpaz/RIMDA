function [pnti,current,fwd] = rimda(M, hs, gain, pnt)
% This function runs RIMDA source localization, finding dipoles for a
% measured field M, a vector of 248 values in Tesla, the number of channels
% for Magnes WH 3600 system.
%   Input:
% M is the 248 by 1 vector representing the measured field
% hs is a N by 3 (xyz) containing the digitized headshape coordinates.
% gain is the channels by locations gain matrix
% pnt is points inside the head arranged in 3 spheres. points outside the
% head are excluded.
%   Output:
% current is the estimated current at the chosen dipoles.
% pnti is the index of the chosen points in the spheric grid pnt.
% fwd is the reconstructed field based on these sources.

% Yuval Harpaz, Oshrit Arviv and Mordehay Medvedovsky
%%
% number of iterations
N = 10000;
%% run iterations
% select 10 locations in each run of the loop and check their possible
% contribution to the measured field M.
Npnt = length(pnt);
% prepare an empty vector Pow to store the results of all the permutations
Pow = zeros(2*Npnt,size(M,2));
% loop fo N times
for permi = 1:N
    % choose 10 locations
    Ran = [];
    [~,ran] = sort(rand(1,Npnt));
    selected = ran(1:10);
    Ran = [Ran;selected];
    srcPerm = false(1,Npnt);
    srcPerm(Ran) = true;
    % use the gain matrix of the chosen 10 points only
    Gain = gain(:,[srcPerm,srcPerm]);
    source = Gain\M; % left divide (similar to pinv)
    recon = Gain*source; % reconstruct the field from hypothesised sources
    R = corr(recon(:),M(:)).^100; % compute correlation between measured and reconstructed field
    % use power of 100 to supress solutions
    % that have less than excellent fit
    pow = zeros(size(Pow,1),size(M,2));
    pow([srcPerm,srcPerm],:) = source.*R; % weight the sources by r^100
    Pow = Pow+pow;    % sum the result with previous ones
end

%% Choose sources
% getCurrent evaluates which sources are 'local maxima' within a radius of
% 30mm, and pass threshold of at least 0.3 strength compared to the strongest source.
% it returns current for sources that are strong maxima. pnti is the index
% of the chosen points. fwd is the reconstructed field based on these
% sources.
[current,~,pnti,fwd] = getCurrent(Pow,pnt,M,gain,30,0.3);
% plot headshape and selected solutions
figure;
plot3(hs(:,1),hs(:,2),hs(:,3),'.k')
hold on
plot3(pnt(pnti,1),pnt(pnti,2),pnt(pnti,3),'o','MarkerSize',9,'MarkerFaceColor','r')
axis equal


