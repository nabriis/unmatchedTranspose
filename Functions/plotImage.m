function [] = plotImage(x,radius,double)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    
    N = length(x);
    N = sqrt(N);
    
    x = reshape(x,N,N);
    
    imagesc(x);
    axis image
    axis off
    colormap gray
    colorbar
    set(gca,'fontsize',16)
    %keyboard
    %set(gcf,'fontsize',16)
    
else
    if nargin<3
        if radius >0
            N = sqrt(length(x));
            xx = reshape(x,N,N);
            %keyboard
            if radius <= N/5; iidx = floor(N/2-radius*2):ceil(N/2+radius*2); end;
            if radius > N/5; iidx = floor(N/2-radius*1.2):ceil(N/2+radius*1.2); end;
            %keyboard
            %xx = xx(iidx,iidx);
            imagesc(xx);
            %keyboard
            colormap gray;
            title('Zoom near ROI')
            axis image
            colorbar
            th = 0:pi/50:2*pi;
            xunit = radius * cos(th) + N/2;
            yunit = radius * sin(th) + N/2;
            hold on; plot(xunit, yunit,'r--','linewidth',1.5);
            hold off;
        end
    else
        N = length(x);
        N = sqrt(N);
        x = reshape(x,N,N);
        imagesc(x);
        axis image
        colormap gray
        colorbar
        
        %hold on
        
        axes('units','normalized','position',[0.01 0.66 0.33 0.33])
        xx = x;
        iidx = floor(N/2-radius):ceil(N/2+radius);
        xx = xx(iidx,iidx);
        imagesc(xx)
        axis off
        %zoom(2)
        
    end
end
end