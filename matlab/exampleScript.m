
cd data
load avgFilt avg1_handL
MHL=avg1_handL.avg(:,138);
load avgFilt avg1_handR
MHR=avg1_handR.avg(:,138);
load avgFilt avg1_footL
MF=avg1_footL.avg(:,180);

M=MF+MHR+MHL;
%rimda(M)

[pnti,current,fwd]=rimda(M);