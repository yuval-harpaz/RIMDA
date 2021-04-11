
% get topography of left hand stimulation field
load avg avg1_handL avg1_handR avg1_footL
MHL = avg1_handL.avg(:,138);
% right hand
MHR = avg1_handR.avg(:,138);
% left foot
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