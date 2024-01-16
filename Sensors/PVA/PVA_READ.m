function [Z_PVA,newPVA,flagPVA] = PVA_READ(PVA_data,mPVA,param)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% PVAmodel = param.PVA.PVAmodel;

PVApos = PVA_data{1}(:,mPVA);
PVAvel = PVA_data{2}(:,mPVA);
PVAatt = PVA_data{3}(:,mPVA);

% switch PVAmodel
%     case 1     
%         Z_PVA = PVApos;
%     case 2
%         Z_PVA = PVAvel;
%     case 3
%         z_PVA = z_PVAatt;
%     case 4
%         z_PVA = [z_PVApos; z_PVAvel];
%     case 5
%         z_PVA = [z_PVApos; z_PVAatt];
%     case 6
%         H_PVA = [H_PVAvel; H_PVAatt];
%         z_PVA = [z_PVAvel; z_PVAatt];
%         R_PVA = blkdiag(R_PVAvel,R_PVAatt);
%     case 7
%         H_PVA = [H_PVApos; H_PVAvel; H_PVAatt];
%         z_PVA = [z_PVApos; z_PVAvel; z_PVAatt];
%         R_PVA = blkdiag(R_PVApos,R_PVAvel,R_PVAatt);
%     otherwise
%         error('Select a case')
% end

Z_PVA = {PVApos,PVAvel,PVAatt};
newPVA = true;
flagPVA = true;%mPVA<=300;

end