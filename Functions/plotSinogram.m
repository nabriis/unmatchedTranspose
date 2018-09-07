function [] = plotSinogram(b,Ntheta)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if length(Ntheta) > 1
    Ntheta = length(Ntheta);
end
p = length(b)/Ntheta;

b = reshape(b,p,Ntheta);

imagesc(b);
%axis equal
axis image
%axis off
colormap gray
colorbar
set(gca,'fontsize',16)
xlabel({'Projection angle'},'interpreter','latex','fontsize',18)
ylabel({'Detector pixel'},'interpreter','latex','fontsize',18)
set(gca,'XTick',[]);
set(gca,'YTick',[]);
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);


end

