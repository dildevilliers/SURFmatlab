spc = [5 1 0.5 0.2 0.1];

for ii = 1:length(spc)
    th = 0:deg2rad(spc(ii)):pi;
    ph = 0:deg2rad(spc(ii)):2*pi;
    [TH, PH] = meshgrid(th,ph);
    R = zeros(size(TH));
    [Fsm0n,Fsmmn,Fsmpn] = FsmnFast(R,TH,PH,20,20,1);
    savestr = strrep(num2str(spc(ii)),'.','p');
    save(['SWEpreBuiltF_',savestr,'_deg_spacing'],'Fsm0n','Fsmmn','Fsmpn')

end