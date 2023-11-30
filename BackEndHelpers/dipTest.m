function [p,dip,xl,xu]=dipTest(x,delta_x)
% Hartigan's dip test of unimodality
% [p, dip, xl, xu] = dipTest (x, delta_x)
% x: vector of observations
% delta_x: optional argument defining the spacing of a discrete distribution
%          (missing delta_x or delta_x=0 means continuous distribution)
% p: p value for rejecting the null hypothesis that the distribution is unimodal
%    Be aware that this p value is not based on an exact calculation but on a look-up table.
%    This table has been constructed using simulations and interpolation is involved
%    when reading it out.
% dip: the dip statistic
% xl: the lower end of the modal interval
% xu: the upper end of the modal interval
%
% J. Ditterich, 12/01

% This is an implementation of the dip test of unimodality published in
% J. A. Hartigan, P. M. Hartigan: The Dip Test of Unimodality. Annals of Statistics,
% 13 (1985) 70-84. It is based on a Fortran algorithm published in
% P. M. Hartigan: Algorithm AS 217: Computation of the Dip Statistic to Test for
% Unimodality. Applied Statistics, 34 (1985) 320-325. It also takes a correction
% into account, which was published in C. J. Sommer, J. N. McNamara: Power considerations
% for the dip test of unimodality using mixtures of normal and uniform distributions.
% American Statistical Association: Proceedings of the Statistical Computing Section,
% 1987, 186-191.
% I have added some extensions:
% 1) Extended table for being on the conservative side even for larger sample sizes
% 2) Extension for handling discrete distributions

% Version history:
% 11/29/01 JD wrote it
% 12/05/01 handling of discrete distributions added

% Check arguments
if (size(x,1)~=1)&(size(x,2)~=1) % not a vector?
    error('DIP_TEST: x must be a vector!');
end;
if length(x)<15 % vector too short?
    error('DIP_TEST: x is too short!');
end;
if nargin<2
    delta_x=0; % default: continuous distribution
end;
if delta_x<0 % negative spacing?
    error('DIP_TEST: delta_x must not be negative!');
end;
% Check plausibility of given spacing
if delta_x>0
    temp=sort(x);
    temp=diff(temp);
    temp=temp(find(temp~=0));
    temp=unique(temp); % get all non-zero spacings
    temp=temp/delta_x;
    
    if ~isempty(find(temp~=round(temp))) % Are there spacings which are not multiples of the given delta_x?
       error('DIP_TEST: The given delta_x is incompatible with x!');
    end;
end;
% Make continuous values from discrete values
if delta_x>0
    x=x+(rand(size(x,1),size(x,2))-.5)*delta_x;
end;
% Sort the vector
x=sort(x);
n=length(x);
% Further checks
if x(1)==x(n) % all observations identical?
    error('DIP_TEST: All observations must not be identical!');
end;
% Computation
low=1; % lower index
high=n; % upper index
dip=1/n;
xl=x(low);
xu=x(high);

% Establish the indices over which combination is necessary for the convex minorant fit
mn=1;

for j=2:n
    mn(j)=j-1;
    
    while 1
        mnj=mn(j);
        mnmnj=mn(mnj);
        a=mnj-mnmnj;
        b=j-mnj;
        
        if (mnj==1)|(((x(j)-x(mnj))*a)<((x(mnj)-x(mnmnj))*b))
            break;
        end;
        
        mn(j)=mnmnj;
    end;
end;

% Establish the indices over which combination is necessary for the concave majorant fit
clear mj;
mj(n)=n;

for jk=1:n-1
    k=n-jk;
    mj(k)=k+1;
    
    while 1
        mjk=mj(k);
        mjmjk=mj(mjk);
        a=mjk-mjmjk;
        b=k-mjk;
        
        if (mjk==n)|(((x(k)-x(mjk))*a)<((x(mjk)-x(mjmjk))*b))
            break;
        end;
        
       mj(k)=mjmjk;
    end;
end;
% Start the cycling
% Collect the change points for the GCM from high to low
clear gcm lcm;
while 1 % line 40 in the Fortran algorithm (big loop)
    ic=1;
    gcm(1)=high;
    
    while 1
        igcm1=gcm(ic);
        ic=ic+1;
        gcm(ic)=mn(igcm1);
        
        if gcm(ic)<=low
            break;
        end;
    end;
    
    icx=ic;
    
    % Collect the change points for the LCM from low to high
    ic=1;
    lcm(1)=low;
    
    while 1
        lcm1=lcm(ic);
        ic=ic+1;
        lcm(ic)=mj(lcm1);
        
        if lcm(ic)>=high
            break;
        end;
    end;
    
    icv=ic;
    ig=icx;
    ih=icv;
    
    % Find the largest distance greater than "dip" between the GCM and the LCM from low to high
    ix=icx-1;
    iv=2;
    d=0;
    
    if (icx~=2)|(icv~=2) % lines 50 - 60 in the Fortran algorithm
        while 1 % line 50 in the Fortran algorithm
            igcmx=gcm(ix);
            lcmiv=lcm(iv);
            
            if igcmx>lcmiv
                % If the next point of either the GCM or LCM is from the GCM then calculate distance here
                lcmiv=lcm(iv); % line 55 in the Fortran algorithm
                igcm=gcm(ix);
                igcm1=gcm(ix+1);
                a=lcmiv-igcm1+1;
                b=igcm-igcm1;
                dx=a/n-((x(lcmiv)-x(igcm1))*b)/(n*(x(igcm)-x(igcm1)));
                iv=iv+1;
                
                if dx>=d
                    d=dx;
                    ig=ix+1;
                    ih=iv-1;
                end;
            else
                % If the next point of either the GCM or LCM is from the LCM then calculate distance here
                lcmiv1=lcm(iv-1);
                a=lcmiv-lcmiv1;
                b=igcmx-lcmiv1-1;
                dx=((x(igcmx)-x(lcmiv1))*a)/(n*(x(lcmiv)-x(lcmiv1)))-b/n;
                ix=ix-1;
                
                if dx>=d
                    d=dx;
                    ig=ix+1;
                    ih=iv;
                end;
            end;
            
            if ix<1 % line 60 in the Fortran algorithm
                ix=1;
            end;
            
            if iv>icv
                iv=icv;
            end;
            
            if gcm(ix)==lcm(iv)
                break; % leave while loop
            end;
        end;
    else
        d=1/n;
    end;
    
    % line 65 in the Fortran algorithm
    if d<dip % Are we done?
        dip=dip/2;
        xl=x(low);
        xu=x(high);
        break; % leave the main loop
    end;
    
    % Calculate the dips for the current low and high
    % The dip for the convex minorant
    dl=0;
    
    if ig~=icx
        icxa=icx-1;
        
        for j=ig:icxa
            temp=1/n;
            jb=gcm(j+1);
            je=gcm(j);
            
            if ((je-jb)>1)&(x(je)~=x(jb))
                a=je-jb;
                const=a/(n*(x(je)-x(jb)));
                
                for jr=jb:je
                    b=jr-jb+1;
                    t=b/n-(x(jr)-x(jb))*const;
                    
                    if t>temp
                        temp=t;
                    end;
                end;
            end;
            
            if dl<temp
                dl=temp;
            end;
        end;
    end;
    
    % The dip for the concave majorant
    du=0;
    
    if ih~=icv
        icva=icv-1;
        
        for k=ih:icva
            temp=1/n;
            kb=lcm(k);
            ke=lcm(k+1);
            
            if ((ke-kb)>1)&(x(ke)~=x(kb))
                a=ke-kb;
                const=a/(n*(x(ke)-x(kb)));
                
                for kr=kb:ke
                    b=kr-kb-1;
                    t=(x(kr)-x(kb))*const-b/n;
                    
                    if t>temp
                        temp=t;
                    end;
                end;
            end;
            
            if du<temp
                du=temp;
            end;
        end;
    end;
    
    % Determine the current maximum
    dipnew=dl;
    
    if du>dl
        dipnew=du;
    end;
    
    if dip<dipnew
        dip=dipnew;
    end;
    
    low=gcm(ig);
    high=lcm(ih);
end; % big loop

% Calculate the p value (based on the table given in the publications)
% I have extended the table for sample sizes of 500 and 1000.
nvec=[15 20 30 50 100 200 500 1000];
dcell=cell(1,8);
dcell{1}=[.0544 .0606 .0641 .0836 .1097 .1179 .1365 .1424 .1538];
dcell{2}=[.0474 .0529 .0569 .0735 .0970 .1047 .1209 .1262 .1382];
dcell{3}=[.0395 .0442 .0473 .0617 .0815 .0884 .1012 .1061 .1177];
dcell{4}=[.0312 .0352 .0378 .0489 .0645 .0702 .0804 .0842 .0926];
dcell{5}=[.0228 .0256 .0274 .0355 .0471 .0510 .0586 .0619 .0987];
dcell{6}=[.0165 .0185 .0197 .0255 .0341 .0370 .0429 .0449 .0496];
dcell{7}=[.0107 .0119 .0127 .0164 .0218 .0237 .0270 .0287 .0311];
dcell{8}=[.0075 .0085 .0091 .0117 .0155 .0168 .0197 .0206 .0233];
p_vec=[.99 .95 .9 .5 .1 .05 .01 .005 .001];

% find nearest n
n_diff=nvec-n;
n_ind=find(n_diff==0);

if isempty(n_ind) % n not tabulated?
    n_ind=find(n_diff>0);
    
    if isempty(n_ind) % Is there no larger n value in the table?
       n_ind=length(nvec); % choose the largest one
    end;
end;
n_ind=n_ind(1); % Choose (1) the exact n if it is in the table;
                %        (2) a value of n which is tabulated, larger than the real n, and has minimum distance;
                %            by choosing a larger n the test will be on the conservative side;
                %        (3) the largest n in the table if the real n is larger than this value.
n_comp=nvec(n_ind); % This is the n for the test.
d_test=dip*sqrt(n/n_comp); % This is the dip value for the test (interpolation based on sqrt(n)*dip as suggested in the paper).
% get the p value
if d_test<min(dcell{n_ind}) % out of range?
    p=p_vec(1); % return the maximum p value in the table
    return;
end;
if d_test>max(dcell{n_ind}) % out of range?
    p=p_vec(length(p_vec)); % return the minimum p value in the table
    return;
end;
% interpolation
p=interp1(dcell{n_ind},p_vec,d_test,'spline');