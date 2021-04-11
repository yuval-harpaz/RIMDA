
% get topography of left hand stimulation field
load avgFilt avg1_handL
MHL = avg1_handL.avg(:,138);
% right hand
load avgFilt avg1_handR
MHR = avg1_handR.avg(:,138);
% left foot
load avgFilt avg1_footL
MF = avg1_footL.avg(:,180);
% merge fields of 3 dipoles
M = MF+MHR+MHL;


% load headshape, gain matrix and shperical grids of points
load hs hs
load gain gain
load pnt pnt
% rimda(M)
rng(2)
[pnti,current,fwd] = rimda(M, hs, gain, pnt);