function [smoothed] = smooth_pad(x,filtsize)
%SMOOTH_PAD Summary of this function goes here
if isempty(x)
    smoothed = x; 
else
    filt_box = (1/filtsize)*ones(1,filtsize); 
    if any(size(x)==1)
    
        x = x(:)'; 
        
        padded = [repmat(x(1),1,filtsize), x, repmat(x(end),1,filtsize)]; 

        filtrd = conv2(padded,filt_box,'same'); 
        filtrd(1:filtsize) = []; 
        filtrd((end-filtsize+1):end) = []; 

        smoothed = filtrd; 
        
    else
        
        padded = [repmat(x(1,:),filtsize,1); x; repmat(x(end,:),filtsize,1)]; 
        
        filtrd = conv2(padded,filt_box','same'); 
        filtrd(1:filtsize,:) = []; 
        filtrd((end-filtsize+1):end,:) = []; 

        smoothed = filtrd; 
    end
        
end

end

