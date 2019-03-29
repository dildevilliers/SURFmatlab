function F = ODReq26(Dm,th_U,th_0)
F = Dm./(4.*(tan(-th_U./2) - tan(-th_0./2)));
end