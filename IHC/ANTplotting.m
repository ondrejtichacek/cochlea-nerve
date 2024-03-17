function [  ] = ANTplotting( SynResFile, NerveResFile )
%ANTPLOTTING 

figure;
surf(SynResFile{1}.V, 'LineStyle', 'none')

figure;
surf(SynResFile{1}.k.Hz, 'LineStyle', 'none')
% figure;
% surf(SynRes(SEL).c, 'LineStyle', 'none')
% figure;
% surf(sum(cat(3,NerveRes(SEL).V),3), 'LineStyle', 'none')
% figure;
% surf(NerveRes(SEL).m, 'LineStyle', 'none')
% figure;
% surf(NerveRes(SEL).h, 'LineStyle', 'none')
% figure;
% surf(NerveRes(SEL).N_Na, 'LineStyle', 'none')
% figure;
% surf(NerveRes(SEL).Istim_HH, 'LineStyle', 'none')

V = NerveResFile{1}.V;
for i = 2:length(NerveResFile)
    V = V + NerveResFile{i}.V;
end

N_Na = NerveResFile{1}.N_Na;
for i = 2:length(NerveResFile)
    N_Na = N_Na + NerveResFile{i}.N_Na;
end

figure
colormap(plasma)
imagesc(V);
set(gca,'YDir','reverse');
colorbar

figure
colormap(plasma)
imagesc(N_Na);
set(gca,'YDir','reverse');
colorbar



end

