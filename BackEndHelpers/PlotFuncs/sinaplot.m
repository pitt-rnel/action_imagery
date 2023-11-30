function [xs,ys,violin_xs,violin_ys] = sinaplot(X,do_plot)
if nargin < 2 
    do_plot = true; 
end
% X is a m-by-1 matrix with m samples
% Sebastien  De Landtsheer November 2018
% sebdelandtsheer@gmail.com

[m, p]=size(X);

MaxDatapoints=sum((~isnan(X))); %The maximum number of non-NaN datapoints in any column
jit=(rand(size(X))-0.5)*0.6; %Uniform jitter centered on 0 
xdef=repmat((1:p),m,1); %Background X-axis position
%Figuring out the density of points
X=sort(X);
Dens=zeros(m,p);
k=2; Dens(k,:)=3./((X(k+1,:)-X(k-1,:)));
k=m-1; Dens(k,:)=3./((X(k+1,:)-X(k-1,:)));
for k=3:m-2
    Dens(k,:)=5./((X(k+2,:)-X(k-2,:)));
end
if MaxDatapoints>16
    for Lim=4:floor((min(40,sqrt(MaxDatapoints)))/2)
        for k=(Lim+1):(m-Lim)
            Dens(k,:)=(Lim*2+1)./((X(k+Lim,:)-X(k-Lim,:)));
        end
    end
end
Dens=Dens./max(Dens(:)); %normalizing density
xval=xdef+(jit.*Dens); %new X-axis values
% xdens1=xdef-(Dens./2);
% xdens2=xdef+(Dens./2);
if do_plot
    figure; hold on; 
end
Colors=lines(10); 
for Group=1:p
    if do_plot
        plot(xval(:,Group), X(:,Group), '.', 'Color', Colors(Group,:)), hold on,
        plot([Group-0.3,Group+0.3],[nanmean(X(:,Group)), nanmean(X(:,Group))],'-k' )
    end
    xs = xval(:,Group);
    ys = X(:,Group); 
end

if do_plot
    set(gca, 'XTick', 1:p);
end

[xv,yv] = violin(X(:,Group),0); 
violin_xs = xv{1}; 
violin_ys = yv{1}; 
if do_plot
    plot(violin_xs,violin_ys,'Color',[.5 .5 .5]); 
end