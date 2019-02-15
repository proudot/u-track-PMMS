function [frameInfo imgDenoised] = main283AUTO_standalone(img, iclean)

% NOTE: this function has been modified only in the way it reads input/returns output.
% The rest of the implementation is identical to main283AUTO.m
% This version operates on a
% The functionality is exactly the same as used in the Loerke09 paper.


% Reconstructs a filtered image using iterative filtering from significant
% coefficients ( see Book by Starck et al, p. 66)

%Then calculates a binary image and the product of the binary image and the
%reconstructed image

% It then finds and plots ALL the local maxima on each cluster (at least
% one maximum per cluster)

%If a cluster has more than one maximum, it divides them into primary and
%secondary (The secondaries are listed at the end.)

%It calculates the intensity of each maximum and the total intensity and size (i.e. area) of
%each cluster and then creates a movie.

%%========================================================================
%
%   NOTE: this is Dinah's modified version, which runs under
%   more recent versions of matlab and therefore on the LINUX server
%
%   Last Date of modification: March 14, 2007
%
%=========================================================================

%   INPUT:  iclean (optional) =     input iclean value; if no value is set
%                                   here, or input is empty, the default
%                                   value will be
%                                   iclean=0
%           icut   (optional) =     input icut value; if no value is set
%                                   here, or input is empty, the default
%                                   value will be
%                                   icut=5
%           itest   (optional) =    input icut value; if no value is set
%                                   here, or input is empty, the default
%                                   value will be
%                                   itest=0
%
%  OUTPUT
%
%  Cluster-specific values:  
%
%     xav, yav   : centroid coordinates for the clusters
%     num        : number of clusters
%
%  Local maxima (each cluster can have several)
%     xmax, ymax : coordinates of local maxima 
%     inn        : intensity of local maxima
%     intot      : intensity of the cluster (sum of all pixels in cluster)
%                : this is repeated for each maxima within the cluster (change in future version?)
%     csize      : number of pixels in each cluster (same as above)
%     lxm        : number of local maxima in cluster (same again)
%     labl       : cluster to which local maximum belong (i.e., cluster label) 
%     nmax       : number of local maxima in the image


% EXAMPLE: main283AUTO(1,[],0);
%
% NOTE: the default values can be changed below in line 58-60, if desired,
% but generally, different parameter values can be used by entering them in
% the commando line instead of changing the code

%=========================================================================
%
%       set default values for parameters
%
%=========================================================================
% if iclean is not defined or if its input is [], set iclean to default
if nargin < 2 || isempty(iclean)
    iclean = 0;
end

% icut and itest are set to default, since we never use any other values
icut = 0;
itest = 0;



[ny,nx,nz] = size(img);
imgDenoised = zeros(ny,nx,nz, 'uint8');


% %%%%%%%%%%%%%%%%%%%%%% SET CONTROL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% icut = 5;
% itest = 0;
% %%%%%%%%%%%%%%%%%%%%%%% SET CONTROL PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% %use iclean = 0 when no only minimal cleaning of very small clusters is necessary (SEE line 163 below)!!!!<<<<<<<<<<(NORMAL SITUATION)
% %use iclean = 1 when some extra cleaning is necessary
% %use iclean = 2 when some more extra cleaning is necessary
% %use iclean = -1 when even minimal cleaning is not necessary
%
% %%% use icut > 0 with Dafi's movie's to remove first few (usually 5) rows (usually contain a lot of noise)
%
% % use itest = 1 to test the goodness of the cluster detection for the first frame
% % use itest = 0 to make a regular detection run for the clusters in all frames

% fprintf(' icut = %d   {icut > 0 for Dafi"s movies to remove first "icut" rows} \n', icut)
% fprintf(' iclean = %d  {iclean = 0 is NORMALLY used} \n', iclean)
% fprintf(' itest = %d {use itest = 0 for a regular run over a whole movie} \n', itest)

% kk=0; %kk is necessary (used as a counter below) (Sylvain)

total_frame_num = size(img, 3);


ixmax = total_frame_num;
if itest == 1
    ixmax = 1;
end

for ix = 1:ixmax
    
    r = img(:,:,ix);
    [sy,sx]=size(r);
    %if icut==1
    r = r(1+icut:sy,1:sx);
    [sy,sx]=size(r);
    %end
    r=double(r);
    mxo=max(r(:));
    mno=min(r(:));
    r=r*450/mxo; % to renormalize the maximum intensity to a reasonable value.
    %avo=mean(r(:));
    %sdo=std(r(:));
    
    k = 4;  %We'll calculate details up to order k
    
    rw=wstarck2(r,k,0); % CALCULATES reconstructed image up to order k without kth detail(irecon=0)
    er=r-rw;  % This is the error image
    
    % Now do iterative filtering
    delta=10;
    n=0;
    sig1=std(er(:));
    while delta > 0.002
        
        n = n+1;
        rw = rw + wstarck2(er,k,0);
        
        er=r-rw;
        sig2=std(er(:));
        delta=abs(sig2-sig1)/sig2;
        sig1=sig2;
    end
    wm = wstarck222(rw,k); % calculate the multiscale product of the details
    rw = rw-min(rw(:));
    
    % mnw=min(rw(:)); (Sylvain)
    
    % size of mask for local average= (2*bs+1)by(2*bs+1)
    [av,sg]=wlocav(rw, 4); % calculate local average (av) and local standard deviation (sg)
    avt = mean(rw(:)); %+0.5*ttsg;
    
    %create binary image
    
    rwb = zeros(sy, sx);
    rwb((rw >= av+0.5*sg) & (rw .* wm >=avt)) = 1;
    
    rwb=bwmorph(rwb,'clean'); %to get rid of isolated pixels (noise)
    rwb=bwmorph(rwb,'fill'); %to fill up empty internal pixels
    rwb=bwmorph(rwb,'thicken'); %to make larger clusters because of the harsh cutoff criteria above
    if iclean >=0
        rwb=bwmorph(rwb,'spur'); % to remove single pixels 8-attached to clusters
        rwb=bwmorph(rwb,'spur');
        rwb=bwmorph(rwb,'clean');% to remove any remaining isolated pixels after applying spur twice
    end
    if iclean > 0 %Value of iclean is set in line 56 above
        rwb=bwmorph(rwb,'erode');%extra cleaning of small spots
        
        if iclean >=2
            rwb=bwmorph(rwb,'spur');%extra cleaning of small spots % for extraextra cleaning
        end
        if iclean >=1
            
            rwb=bwmorph(rwb,'clean');%extra cleaning of small spots
        end
        rwb=bwmorph(rwb,'thicken');%extra cleaning of small spots
    end

    [Lbr,num] = bwlabel(rwb,8); %Label the clusters. Better to use label 8 than label 4
    
    rw=rw*(mxo-mno)/max(rw(:)); % normalize maximum intensity to original (maximum intensity-background intensity).
    rw=rw.*rwb; %Outside the clusters the intensity is set to zero exactly
    
    
    lm = locmax2d(rw, [9 9]); % find the location of the maxima of the clusters
    [ymx, xmx] = find(lm>0); % coordinates of the local maxima
    nmax = length(xmx);  % calculate  number of maxima determined by locmax2d
    
    snum=0; %initialize number of secondary maxima (if there are two or more maxima in a cluster found by locmax2d)
    lxm=zeros(num,1);  %initialize the number of maxima in each cluster
    intot=zeros(num,1);%initialize the total intensity of each cluster
    csize=zeros(num,1);%initialize the size of each cluster
    xmax=zeros(num,1);% initialize the coordinates of the maxima
    ymax=zeros(num,1);
    yav=zeros(num,1);%initialize the coordinates of the center of intensity
    xav=zeros(num,1);
    labl=zeros(num,1); %NEWHRJ
    nn=zeros(num, 1); % (Sylvain)
    ymax2=zeros(num, 1); % (Sylvain)
    xmax2=zeros(num, 1); % (Sylvain)
    lxm2=zeros(num, 1); % (Sylvain)
    intot2=zeros(num, 1); % (Sylvain)
    csize2=zeros(num, 1); % (Sylvain)
    labl2=zeros(num, 1); % (Sylvain)
    for i=1:num
        labl(i)=i;%NEWHRJ
        [yi,xi] = find(Lbr==i);
        %kk=kk+1; (Sylvain)
        lni=length(yi);
        nn(i)=lni; %nn(kk)=lni; (Sylvain)
        [ym,xm]= find((Lbr==i) & (lm>0));
        
        lxm(i)=length(xm);

        % pnum=pnum+1; (Sylvain)
        if lxm(i)==1
            
            ymax(i)=ym; % ymax(pnum)=ym; (Sylvain)
            xmax(i)=xm; % xmax(pnum)=xm; (Sylvain)
            for j=1:lni
                inloc=rw(yi(j),xi(j));
                yav(i)=yav(i)+yi(j)*inloc;
                xav(i)=xav(i)+xi(j)*inloc;
                intot(i)=intot(i)+inloc;
                csize(i)=csize(i)+1;
            end
            xav(i)=xav(i)/intot(i);
            yav(i)=yav(i)/intot(i);
        elseif lxm(i)==0 %if the cluster's maximum has not been located by locmax, search for it locally
            nmax=nmax+1;
            
            tt=0;
            lxm(i)=1;
            for j=1:lni
                inloc=rw(yi(j),xi(j));
                yav(i)=yav(i)+yi(j)*inloc;
                xav(i)=xav(i)+xi(j)*inloc;
                intot(i)=intot(i)+inloc;
                %                 intot(1)
                csize(i)=csize(i)+1;
                if inloc>tt
                    tt=inloc;
                    ymx(nmax)=yi(j);
                    xmx(nmax)=xi(j);
                    ymax(i)=yi(j); % ymax(pnum)=yi(j); (Sylvain)
                    xmax(i)=xi(j); % xmax(pnum)=xi(j); (SYlvain)
                end
                %fprintf([' i = ',num2str(i),' j=  ',num2str(j),'\n'])
                %ymax(pnum)=ymx(nmax);
                %xmax(pnum)=xmx(nmax);
            end
            xav(i)=xav(i)/intot(i);
            yav(i)=yav(i)/intot(i);
        else  % Take care of the case when there are more than one maximum in each cluster
            
            for j=1:lni
                inloc=rw(yi(j),xi(j));
                yav(i)=yav(i)+yi(j)*inloc;
                xav(i)=xav(i)+xi(j)*inloc;
                intot(i)=intot(i)+inloc;
                csize(i)=csize(i)+1;
            end
            xav(i)=xav(i)/intot(i);
            yav(i)=yav(i)/intot(i);
            tt1=rw(ym(1),xm(1));
            ymax(i)=ym(1); % ymax(pnum)=ym(1); (Sylvain)
            xmax(i)=xm(1); % xmax(pnum)=xm(1); (Sylvain)
            for jmx=2:lxm(i)
                tt2=rw(ym(jmx),xm(jmx));
                snum=snum+1;
                if tt2>tt1
                    ymax2(snum)=ymax(i); % ymax2(snum)=ymax(pnum); (Sylvain) %The weaker maxima  are designated secondary maxima
                    xmax2(snum)=xmax(i); % xmax2(snum)=xmax(pnum); (Sylvain)
                    ymax(i)=ym(jmx); % ymax(pnum)=ym(jmx); (Sylvain) %The maximum with the highest intensity in a cluster is designated the primary max
                    xmax(i)=xm(jmx); % xmax(pnum)=xm(jmx); (Sylvain)
                    lxm2(snum)=lxm(i);
                    intot2(snum)=intot(i);
                    csize2(snum)=csize(i);
                    labl2(snum)=i;%NEWHRJ
                    tt1=tt2;
                else
                    ymax2(snum)=ym(jmx);
                    xmax2(snum)=xm(jmx);
                    lxm2(snum)=lxm(i);
                    intot2(snum)=intot(i);
                    csize2(snum)=csize(i);
                    labl2(snum)=i;%NEWHRJ
                end
            end
        end
%         return;
    end
    % Add secondary maxima to the end of the primary maxima list
    if snum > 0
        % REPLACE THIS... (Sylvain)
        %         for  imx=1:snum
        %             pnum=pnum+1;
        %             ymax(pnum)=ymax2(imx);
        %             xmax(pnum)=xmax2(imx);
        %             lxm(pnum)=lxm2(imx);
        %             intot(pnum)=intot2(imx);
        %             csize(pnum)=csize2(imx);
        %             labl(pnum)=labl2(imx);%NEWHRJ
        %         end
        % ...BY THIS:
        ymax = vertcat(ymax, ymax2(1:snum));
        xmax = vertcat(xmax, xmax2(1:snum));
        lxm = vertcat(lxm, lxm2(1:snum));
        intot = vertcat(intot, intot2(1:snum));
        csize = vertcat(csize, csize2(1:snum));
        labl = vertcat(labl, labl2(1:snum));
    end
    nmax=length(xmax);
        
    inn = zeros(nmax,1);
    for ii = 1:nmax
        inn(ii) = rw(ymax(ii),xmax(ii)); % problem with shift? xmax, ymax start at 0?     
    end
    
    frameInfo(ix).ymax = ymax;
    frameInfo(ix).xmax = xmax;
    frameInfo(ix).inn = inn;
    frameInfo(ix).yav = yav;
    frameInfo(ix).xav = xav;
    frameInfo(ix).intot = intot;
    frameInfo(ix).csize = csize;
    frameInfo(ix).lxm = lxm;
    frameInfo(ix).labl = labl;
    frameInfo(ix).num = num;
    frameInfo(ix).nmax = nmax;
   
    % eliminate the frequent error messages by rescaling to uint8
    rw8 = uint8(round(255*(rw/max(rw(:)))));
    imgDenoised(:,:,ix) = rw8;
end
