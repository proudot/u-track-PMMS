function [pValues,distributions,limits]=chi2correlation(XY, cutoff, limits, expect, verbose)
%CHI2CORRELATION calculates the correlation between two variables
%
% SYNOPSIS: [pValues]=chi2correlation(XY)
%
% INPUT XY:           n-by-2 array of the two variables to correlate
%                     Alternatively, specify handles to axes with the data
%                     If there are less than 50 data points, the code will
%                     return a warning.
%       cutoff (opt): % of data that should be removed in each dimension.
%                     Cutoff of 10 will remove the smallest and largest 5%
%                     of X and Y, respectively
%                     Alternatively, specify a 4-element vector with
%                     [xmin,xmax,ymin,ymax]
%                     Default: 0
%       limits (opt): vector with [xmin,xmax,ymin,ymax] as the limits for
%                     data. If not specified, the minimum and maximum x and
%                     y values, respectively, will be used as limits.
%       expect (opt): 4-by-7 array with expected distribution of the data
%                     (see also description of output). Specify NaN if you
%                     want to use the default distributions.
%                     Default: all NaN
%       verbose (opt): Whether to plot a figure. Default: 1
%
% OUTPUT pValues:     1-by-7 vector of p-values corresponding to the
%                     following tests (all chi2-tests)
%                     1) independence of X and Y, using the four quadrants
%                     2) uniform distribution of the data among four
%                        equally spaced horizontal bands
%                     3) uniform distribution of the data among four
%                        equally spaced vertical bands
%                     4,5) independence of X and Y using four equally
%                          spaced diagonal bands ordered from the bottom
%                          left corner to the top right corner
%                          (significance indicates that there may be a
%                          cutoff of the form x+y = c).
%                          4: the slices are handled as if they were of the
%                             same size (tests expected counts)
%                          5: the corner slices are increased by a factor
%                             3 to eliminate size effects (tests expected
%                             densities)
%                     6,7) independence of X and Y using four equally
%                          spaced diagonal bands ordered from the top
%                          left corner to the bottom right corner
%                          (significance indicates that there may be a
%                          cutoff of the form x-y = c)
%                          6: the slices are handled as if they were of the
%                             same size (tests expected counts)
%                          7: the corner slices are increased by a factor
%                             3 to eliminate size effects (tests expected
%                             densities)
%
%       distributions : 4-by-7 matrix of relative distributions of the data
%                       among the four groups that are used in the test
%       limits        : [xmin, xmax, ymin, ymax] as in input of the same
%                       name
%
% REMARKS The lines on top of the bar plot indicate the expected
%         distributions.
%
%         If you know theoretical limits to your distribution, you should
%         use them. Otherwise, the diagonal tests may not show the proper
%         results
%         Example:
%               xy = rand(100,2);
%               xy(sum(xy,2)>1,:)=[];
%               chi2correlation(xy) % data-derived limits
%               chi2correlation(xy,[],[0,1,0,1]) % correct limits
%               chi2correlation(xy,[],[0,1,0,2]) % way off limits
%
%         There should be at least ~50 data points in your comparison. The
%         more data there is, the higher the significance. As an example,
%         repeat the above test with ~500 data points (xy=rand(1000,2)).
%
%
% created with MATLAB ver.: 7.3.0.267 (R2006b) on Windows_NT
%
% created by: jdorn
% DATE: 19-Oct-2006
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%==========================
%% DEFAULTS
%==========================

% cutoff for removing outliers. Cutoff = 10 means the minimum 5% and the
% maximum 5% of both x and y get removed
def_cutoff = 0;


%===========================
%% TEST INPUT
%===========================

if nargin < 1 || isempty(XY)
    error('CHI2CORRELATION:WRONGINPUT',...
        'chi2correlation requres two nonempty input arguments')
elseif length(XY) == 1 && ishandle(XY)
    ah = XY;
    ph = get(ah,'Children');
    % in case there are multiple lines, take the last kid (=first drawn)
    ph = ph(end);
    if ~strcmpi(get(ph,'Type'),'line')
        error('CHI2CORRELATION:WRONGINPUT',...
            'chi2correlation requires that the axis has a lineseries as children')
    end
    X = get(ph,'xData')';
    Y = get(ph,'yData')';
    % read from axes handle with rectangle
    set(ah,'xColor','r','yColor','r')
    def_cutoff = getrect(ah);
    set(ah,'xColor','k','yColor','k')
    % make cutoff [xmin, xmax ymin ymax]
    def_cutoff = def_cutoff([1,3,2,4]);
    % assign default cutoff to potentially allow override
    def_cutoff([2,4]) = abs(def_cutoff([2,4]))+def_cutoff([1,3]);
else
    XY = returnRightVector(XY, 2);
    X = XY(:,1);
    Y = XY(:,2);
end

if ~all(size(X)==size(Y))
    error('CHI2CORRELATION:WRONGINPUT',...
        'x and y have to be of the same size for chi2correlation')
end

if nargin < 2 || isempty(cutoff)
    cutoff = def_cutoff;
else
    % no testing yet
end

if nargin < 3 || isempty(limits)
    limits = [];
end

if nargin < 4 || isempty(expect)
    expect = repmat(NaN,4,7);
end

if nargin < 5 || isempty(verbose)
    verbose = true;
end

%===========================


%===========================
%% PREPARE THE DATA
%===========================

% remove NaNs
goodRows = ~isnan(X) & ~isnan(Y);
X = X(goodRows);
Y = Y(goodRows);

% remove data - either via %cutoff, or via rectangle
if length(cutoff) == 1

    % remove n% of the data
    xCut = prctile(X,[cutoff/2,100-cutoff/2]);
    yCut = prctile(Y,[cutoff/2,100-cutoff/2]);

else
    % explicit cutoff
    xCut = cutoff(1:2);
    yCut = cutoff(3:4);
end

% if there are outside limits: check that cutoff is within
if ~isempty(limits)
    xCut(1) = min(limits(1),xCut(1));
    xCut(2) = max(limits(2),xCut(2));
    yCut(1) = min(limits(3),yCut(1));
    yCut(2) = max(limits(4),yCut(2));
end

% remove first the Y values that go with the x-cutoff and vice versa
Y(X<xCut(1)) = [];
X(X<xCut(1)) = [];
Y(X>xCut(2)) = [];
X(X>xCut(2)) = [];
X(Y<yCut(1)) = [];
Y(Y<yCut(1)) = [];
X(Y>yCut(2)) = [];
Y(Y>yCut(2)) = [];


% scale to 0...1
if isempty(limits)
    
    xmin = min(X);
    X = X - xmin;
    xmax = max(X);
    X = X / xmax;

    ymin = min(Y);
    Y = Y - ymin;
    ymax = max(Y);
    Y = Y / ymax;

else
    xmin = limits(1);
    X = X - xmin;
    xmax = limits(2)-limits(1);
    X = X / xmax;

    ymin = limits(3);
    Y = Y - ymin;
    ymax = limits(4)-limits(3);
    Y = Y / ymax;
end

%===========================


%===========================
%% CALCULATE STATISTICS
%===========================

% general: Values that lie exactly on a line (such as X==0.5) count half on
% each side

nTot = length(X);
% warn if not enough points
if nTot < 50
    warning('CHI2CORRELATION:NOTENOUGHDATA',...
        'chi2correlation will not be particularly sensitive with less than 50 data points (%i)',nTot);
end

% chi-square : split data into four quarters, calculate expected probabilites,
% calculate the chi2-score

% ratios in the four quadrants
Q = zeros(2);
Q(1) = (sum(Y>0.5 & X<0.5) + 0.5*sum(Y==0.5 & X<0.5) + 0.5*sum(Y>0.5 & X==0.5))/nTot;
Q(2) = (sum(Y<0.5 & X<0.5) + 0.5*sum(Y==0.5 & X<0.5) + 0.5*sum(Y<0.5 & X==0.5))/nTot;
Q(3) = (sum(Y>0.5 & X>0.5) + 0.5*sum(Y==0.5 & X>0.5) + 0.5*sum(Y>0.5 & X==0.5))/nTot;
Q(4) = (sum(Y<0.5 & X>0.5) + 0.5*sum(Y==0.5 & X>0.5) + 0.5*sum(Y<0.5 & X==0.5))/nTot;

% top/bottom, left/right
tb = zeros(2);
lr = zeros(2);
tb(1,:) = Q(1) + Q(3);
tb(2,:) = Q(2) + Q(4);
lr(:,1) = Q(1) + Q(2);
lr(:,2) = Q(3) + Q(4);

% expected values
if any(isnan(expect(:,1)))
    EQ = tb .* lr;
else
    % use outside expectation
    EQ = expect(:,1);
end

% prevent zeros
EQ = max(EQ,realmin);

% chi2-score, p-value
chi2chi = sum(sum((Q(:)-EQ(:)).^2./EQ(:)))*nTot;
pChi = 1-chi2cdf(chi2chi,1);


% vertical stripes (cols)
V = zeros(4,1);
V(1) = (sum(X < 0.25) + 0.5*sum(X == 0.25))/nTot;
V(2) = (sum(X > 0.25 & X < 0.5) + 0.5*sum(X == 0.25) + 0.5*sum(X == 0.5))/nTot;
V(3) = (sum(X > 0.5  & X < 0.75) + 0.5*sum(X == 0.5) + 0.5*sum(X == 0.75))/nTot;
V(4) = (sum(X > 0.75) + 0.5*sum(X == 0.75))/nTot;

% expectation: homogeneous distribution
if any(isnan(expect(:,3)))
    EV = repmat(0.25,4,1);
else
    % use outside expectation
    EV = expect(:,3);
end

% prevent zeros
EV = max(EV,realmin);

% chi2-score, p-value
chi2vert = sum((V-EV).^2./EV)*nTot;
pVert = 1-chi2cdf(chi2vert,3);


% horizontal stripes (rows)
% !! we go top-to-bottom

H = zeros(4,1);
H(4) = (sum(Y < 0.25) + 0.5*sum(Y == 0.25))/nTot;
H(3) = (sum(Y > 0.25 & Y < 0.5) + 0.5*sum(Y == 0.25) + 0.5*sum(Y == 0.5))/nTot;
H(2) = (sum(Y > 0.5 & Y < 0.75) + 0.5*sum(Y == 0.5) + 0.5*sum(Y == 0.75))/nTot;
H(1) = (sum(Y > 0.75) + 0.5*sum(Y == 0.75))/nTot;

% expectation: homogeneous distribution
if any(isnan(expect(:,2)))
    EH = repmat(0.25,4,1);
else
    % use outside expectation
    EH = expect(:,2);
end

% prevent zeros
EH = max(EH,realmin);

% chi2-score, p-value
chi2horz = sum((H-EH).^2./EH)*nTot;
pHorz = 1-chi2cdf(chi2horz,3);


% multiply Vert and Horz
VH = repmat(V',4,1) .* repmat(H,1,4);


% diagonal stripes BottomLeft to TopRight
% X+Y == 1 is the line from TopLeft to BottomRight

BLTR = zeros(4,1);
BLTR(1) = (sum(X + Y < 0.5) +           + 0.5*sum(X + Y == 0.5)                      )/nTot;
BLTR(2) = (sum(X + Y > 0.5 & X + Y < 1) + 0.5*sum(X + Y == 0.5) + 0.5*sum(X + Y == 1))/nTot;
BLTR(3) = (sum(X + Y < 1.5 & X + Y > 1) + 0.5*sum(X + Y == 1.5) + 0.5*sum(X + Y == 1))/nTot;
BLTR(4) = (sum(X + Y > 1.5) +             0.5*sum(X + Y == 1.5)                      )/nTot;

% expectation: Calculate from horizontal and vertical stripes under the
% assumption of independent variables
d1 = diag(1,-3);
d2 = diag([0.5,0.5],-2);
d3 = diag([1,1,1],-1);
d4 = diag([0.5,0.5,0.5,0.5],0);
d5 = diag([1,1,1],1);
d6 = diag([0.5,0.5],2);
d7 = diag(1,3);

if any(isnan(expect(:,4)))
    EBLTR = zeros(4,1);
    EBLTR(1) = sum(sum((d1 + d2).*VH));
    EBLTR(2) = sum(sum((d2 + d3 + d4).*VH));
    EBLTR(3) = sum(sum((d4 + d5 + d6).*VH));
    EBLTR(4) = sum(sum((d6 + d7).*VH));
else
    % use outside expectation
    EBLTR = expect(:,4);
end

% prevent zeros
EBLTR = max(EBLTR,realmin);

% chi2-score, p-value
chi2bltr = sum((BLTR-EBLTR).^2./EBLTR)*nTot;
pBLTR = 1-chi2cdf(chi2bltr,3);

% chi2-score, p-value weighted
BLTRw = BLTR.*[3;1;1;3];
BLTRw = BLTRw/sum(BLTRw);
if any(isnan(expect(:,5)))
    EBLTRw = EBLTR.*[3;1;1;3];
    EBLTRw = EBLTRw/sum(EBLTRw);
else
    % use outside expectation
    EBLTRw = expect(:,5);
end
chi2bltrw = sum((BLTRw-EBLTRw).^2./EBLTRw)*nTot;
pBLTRw = 1-chi2cdf(chi2bltrw,7);


% diagonal stripes TopLeft to BottomRight
% X-Y == 0 is the line from BottomLeft to TopRight

TLBR = zeros(4,1);
TLBR(1) = (sum(X - Y < -0.5) +                + 0.5*sum(X - Y == -0.5)                      )/nTot;
TLBR(2) = (sum(X - Y > -0.5 & X - Y < 0) + 0.5*sum(X - Y == -0.5) + 0.5*sum(X - Y == 0))/nTot;
TLBR(3) = (sum(X - Y <  0.5 & X - Y > 0) + 0.5*sum(X - Y ==  0.5) + 0.5*sum(X - Y == 0))/nTot;
TLBR(4) = (sum(X - Y >  0.5) +                  0.5*sum(X - Y ==  0.5)                      )/nTot;

% expectation: Calculate from horizontal and vertical stripes under the
% assumption of independent variables
d1 = d1(end:-1:1,:);
d2 = d2(end:-1:1,:);
d3 = d3(end:-1:1,:);
d4 = d4(end:-1:1,:);
d5 = d5(end:-1:1,:);
d6 = d6(end:-1:1,:);
d7 = d7(end:-1:1,:);

if any(isnan(expect(:,6)))
    ETLBR = zeros(4,1);
    ETLBR(1) = sum(sum((d1 + d2).*VH));
    ETLBR(2) = sum(sum((d2 + d3 + d4).*VH));
    ETLBR(3) = sum(sum((d4 + d5 + d6).*VH));
    ETLBR(4) = sum(sum((d6 + d7).*VH));
else
    % use outside expectation
    ETLBR = expect(:,6);
end

% prevent zeros
ETLBR = max(ETLBR,realmin);


% chi2-score, p-value
chi2tlbr = sum((TLBR-ETLBR).^2./ETLBR)*nTot;
pTLBR = 1-chi2cdf(chi2tlbr,3);

% chi2-score, p-value weighted
TLBRw = TLBR.*[3;1;1;3];
TLBRw = TLBRw/sum(TLBRw);
if any(isnan(expect(:,7)))
    ETLBRw = ETLBR.*[3;1;1;3];
    ETLBRw = ETLBRw/sum(ETLBRw);
else
    % use outside expectation
    ETLBRw = expect(:,7);
end

chi2tlbrw = sum((TLBRw-ETLBRw).^2./ETLBRw)*nTot;
pTLBRw = 1-chi2cdf(chi2tlbrw,3);

%======================



%======================
%% PLOT
%======================
if verbose
    figure
    % plot data
    ah = subplot(2,3,1);
    h=scattercloud(X*xmax+xmin,Y*ymax+ymin,[],[],'b.');
    xlim([xmin,xmin+xmax])
    ylim([ymin,ymin+ymax])
    hold on
    line([(2*xmin+xmax)/2;(2*xmin+xmax)/2;NaN;xmin;xmin+xmax;NaN;xmin;xmin+xmax;NaN;xmin;xmin+xmax],...
        [ymin;ymin+ymax;NaN;(2*ymin+ymax)/2;(2*ymin+ymax)/2;NaN;ymin;ymin+ymax;NaN;ymin+ymax;ymin],'Color','r')
    xlabel('L -> R')
    ylabel('B -> T')
    
    % set point size of scattercloud to 3
    set(h(2),'MarkerSize',1)
    set(get(ah,'Title'),'String',sprintf('N = %i',nTot))


    xx = [0.125,0.375,0.625,0.875]';

    % plot squares
    ah = subplot(2,3,4);
    bh=bar(xx,Q(:),1);
    set(bh,'FaceColor','b','EdgeColor','none');
    hold on,
    bh=bar(xx,EQ(:),1);
    set(bh,'FaceColor','none','EdgeColor','r','LineWidth',1);
    ylim([0,1])
    set(ah,'xTick',xx,'xTickLabel',['Q1';'Q2';'Q3';'Q4']);
    set(get(ah,'Title'),'String',sprintf('pChi=%1.4f',pChi))

    % plot horizontal
    ah = subplot(2,3,2);
    bh=bar(xx,H,1);
    set(bh,'FaceColor','b','EdgeColor','none');
    hold on,
    bh=bar(xx,EH(:),1);
    set(bh,'FaceColor','none','EdgeColor','r','LineWidth',1);
    ylim([0,1])
    set(ah,'xTick',xx,'xTickLabel',['1';'2';'3';'4']);
    set(get(ah,'Title'),'String',sprintf('pHorz (T->B)=%1.4f',pHorz))

    % plot vertical
    ah = subplot(2,3,5);
    bh=bar(xx,V,1);
    set(bh,'FaceColor','b','EdgeColor','none');
    hold on,
    bh=bar(xx,EV(:),1);
    set(bh,'FaceColor','none','EdgeColor','r','LineWidth',1);
    ylim([0,1])
    set(ah,'xTick',xx,'xTickLabel',['1';'2';'3';'4']);
    set(get(ah,'Title'),'String',sprintf('pVert (L->R)=%1.4f',pVert))

    % plot BL-TR

    ah = subplot(2,3,3);
    bh=bar(xx,BLTRw,1);
    set(bh,'FaceColor','b','EdgeColor','none');
    hold on,
    bh=bar(xx,EBLTRw(:),1);
    set(bh,'FaceColor','none','EdgeColor','r','LineWidth',1);
    ylim([0,1])
    set(ah,'xTick',xx,'xTickLabel',['1';'2';'3';'4']);
    set(get(ah,'Title'),'String',sprintf('pBL->TR cts:%1.4f rho:%1.4f',pBLTR,pBLTRw))

    % plot TLBR
    ah = subplot(2,3,6);
    bh=bar(xx,TLBRw,1);
    set(bh,'FaceColor','b','EdgeColor','none');
    hold on,
    bh=bar(xx,ETLBRw(:),1);
    set(bh,'FaceColor','none','EdgeColor','r','LineWidth',1);
    ylim([0,1])
    set(ah,'xTick',xx,'xTickLabel',['1';'2';'3';'4']);
    set(get(ah,'Title'),'String',sprintf('pTL->BR cts:%1.4f rho:%1.4f',pTLBR,pTLBRw))
end
%========================


% assign output
if nargout > 0
    pValues = [pChi,pHorz,pVert,pBLTR,pBLTRw,pTLBR,pTLBRw];
    if nargout > 1
        distributions = [Q(:),H,V,BLTR,BLTRw,TLBR,TLBRw];
        limits = [xmin,xmax+xmin,ymin,ymax+ymin];
    end
end

