function [pnti,current,fwd]=rimda(M)
% This function runs RIMDA source localization, finding dipoles for a
% measured field M, a vector of 248 values in Tesla, the number of channels
% for Magnes WH 3600 system.
% We run the function in subject directory where the headshape can be found
% (hs_file), along with data (c,rfhp0.1Hz) and the config file. for the
% current example we suply averaged and filtered data (using FieldTrip),
% headshape points, gain matrix and coordinates of points inside the head
% in matlab format.

% Output
% current is the current at the chosen dipoles. pnti is the index
% of the chosen points in the spheric grid. fwd is the reconstructed field based on these
% sources.
% PowCur needs fixing, it is for multiple timepoints, not applicable for
% one time point.
%%
% number of iterations
N=10000;
% set a spherical grid of equally distributed points (642). here we rely on
% GridSphere, an external package, although FieldTrip creates such a grid
% by default when no grid is specified for creating leadfield matrices.


load hs % digitized headshape
load gain % gain matrix
load pnt % points inside the head arranged in 3 spheres. points outside the head are excluded
   
%% RIMDA stage
% run permutations: select 10 locations in each run of the loop and check
% their possible contribution to the measured field M.
Npnt=length(pnt);
% prepare an empty vector Pow to store the results of all the permutations
Pow=zeros(2*Npnt,size(M,2));
% loop fo N times (10000)
for permi=1:N
    % choose 10 locations
    Ran=[];
    [~,ran]=sort(rand(1,Npnt));
    selected=ran(1:10);
    Ran=[Ran;selected];
    srcPerm=false(1,Npnt);
    srcPerm(Ran)=true;
    % use the gain matrix of the chosen 10 points only
    Gain=gain(:,[srcPerm,srcPerm]);
    source=Gain\M; % left divide (similar to pinv)
    recon=Gain*source; % reconstruct the field from hypothesised sources
    R=corr(recon(:),M(:)).^100; % compute correlation between measured and reconstructed field
                                % use power of 100 to supress solutions
                                % that have less than excellent fit
    pow=zeros(size(Pow,1),size(M,2)); 
    pow([srcPerm,srcPerm],:)=source.*R; % weight the sources by r^100
    Pow=Pow+pow;    % sum the result with previous ones
    %prog(permi) % feedback on progress
end

%% Choose sources
% this function evaluates which sources are 'local maxima' within a radius of  
% 30mm, and pass threshold of at least 0.3 strength compared to the strongest source.
% it returns current for sources that are strong maxima. pnti is the index
% of the chosen points. fwd is the reconstructed field based on these
% sources.
[current,~,pnti,fwd]=getCurrent(Pow,pnt,M,gain,30,0.3,false);
figure;
plot3(hs(:,1),hs(:,2),hs(:,3),'.k')
hold on
plot3(pnt(pnti,1),pnt(pnti,2),pnt(pnti,3),'o','MarkerSize',9,'MarkerFaceColor','r')
% make some plots and prepare output

