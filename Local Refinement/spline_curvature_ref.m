function Xnew = spline_curvature_ref(X, nI0, Q, ref, Ufac, TEfac);
% function spline_curvature_ref(X, nI0, Q, ref);
%
% INPUT
% X   =  points to spline, not including TE
% nI0 =  number of intervals (e.g. 20)
% Q   =  number of points per interval (geometry order, e.g. 3)
% ref =  refinement index (nI = 2^ref*nI0) ... e.g. 2 or 3
% Ufac  = uniformity factor (1 = normal; > 1 means more uniform distribution of points)
% TEfac = trailing-edge resolution factor (1 = normal; > 1 = high; < 1 = low)
%
% OUTPUT
% Xnew = new points (TE not included)
% 

% include TE point at end twice
[~,I] = max(X(:,1)); % ensures always takes TE
X = X([1:size(X,1),I],:);
X(end,1) = X(end,1)+1e-6; % ensures will always work

% min/max of given points (x-coordinate)
xmin = min(X(:,1));
xmax = max(X(:,1));

% spline given points
PPLE = spline2d(X);

% curvature-based spacing on geom
nfine = 500;
stot = PPLE.breaks(end); 
s = linspace(0,stot,nfine);
xyfine = myppval(PPLE, s);
PPfine = spline2d(xyfine);
sk(1) = 0;
[xq, wq] = quad1;
for i = 1:nfine-1,
  ds = s(i+1)-s(i);
  st  = 0.5*(1+xq)*ds;
  px = PPfine.xcoefs(i,:);
  xs = 3.0*px(1)*st.^2 + 2.0*px(2)*st + px(3);
  xss = 6.0*px(1)*st + 2.0*px(2);
  py = PPfine.ycoefs(i,:);
  ys = 3.0*py(1)*st.^2 + 2.0*py(2)*st + py(3);
  yss = 6.0*py(1)*st + 2.0*py(2);
  skint = 0.01*Ufac+0.5*wq'*sqrt(xss.*xss + yss.*yss)*ds;
  
  % force TE resolution
  xx = (0.5*(xyfine(i,1)+xyfine(i+1,1))-xmin)/(xmax-xmin); % close to 1 means at TE
  skint = skint + TEfac*0.5*exp(-100*(1.0-xx));

  % increment sk
  sk(i+1) = sk(i) + skint;
end

% offset by fraction of average to avoid problems with zero curvature
sk = sk + 2.0*sum(sk)/nfine;

% points on nI0 intervals
PPsk = spline(sk,s);

sk = linspace(PPsk.breaks(1),PPsk.breaks(end), nI0+1);

% corresponding s values
s0 = ppval(PPsk, sk);

% account for Q and ref
nsub = 2^ref*Q;
if (nsub == 1)
  s = s0;
else
  nI = nsub*nI0;
  s = zeros(1, nI+1);
  for i=0:nI0-1,
    for j = 1:nsub,
      r = (j-1.0)/(nsub);
      s(i*nsub+j) = s0(i+1)*(1.0-r) + s0(i+2)*r;
    end
  end
end

% new points
Xnew = myppval(PPLE, s);
Xnew = Xnew(1:end-1, :); % do not include TE

%figure(1); clf; hold on;
%plot(Xnew(:,1), Xnew(:,2), 'b.');
%hold off;


%-----------------------------------------------------------
% splines 2d (x,y) points
function [PP] = spline2d(X);
N = length(X);
S    = zeros(N,1);
Snew = zeros(N,1);

% estimate arc-length parameter
S(1) = 0.0;
for i=2:N, 
  S(i) = S(i-1) + sqrt((X(i,1)-X(i-1,1))^2 + (X(i,2)-X(i-1,2))^2);
end

PPX = spline(S,X(:,1));
PPY = spline(S,X(:,2));

% re-integrate to true arc-length via several passes
[xq, wq] = quad1;
for ipass = 1:10,
  serr = 0;
  Snew(1) = S(1);
  for i = 1:(N-1),
    ds = S(i+1)-S(i);
    st  = 0.5*(1+xq)*ds;
    px = PPX.coefs(i,:);
    xs = 3.0*px(1)*st.^2 + 2.0*px(2)*st + px(3);
    py = PPY.coefs(i,:);
    ys = 3.0*py(1)*st.^2 + 2.0*py(2)*st + py(3);
    sint = 0.5*wq'*sqrt(xs.*xs + ys.*ys)*ds;
    serr = max(serr, abs(sint-ds));
    Snew(i+1) = Snew(i) + sint;
  end
  S = Snew;
  PPX = spline(S,X(:,1));
  PPY = spline(S,X(:,2));
end

PP = two2one(PPX, PPY);


%-----------------------------------------------------------
% splits 2d spline into two
function [PP1, PP2] = split_splines(PP, n);

[PPX, PPY] = one2two(PP);
[PPX1, PPX2] = split_spline(PPX,n);
[PPY1, PPY2] = split_spline(PPY,n);
PP1 = two2one(PPX1, PPY1);
PP2 = two2one(PPX2, PPY2);

%-----------------------------------------------------------
% splits 1d spline into two
function [PP1, PP2] = split_spline(PP, n);

nt = PP.pieces+1;
PP1 = PP;
PP2 = PP;
PP1.pieces = n-1;
PP2.pieces = nt-n;
PP1.breaks = PP.breaks(1:n);
PP2.breaks = PP.breaks(n:nt)-PP.breaks(n);
PP1.coefs  = PP.coefs(1:(n-1),:);
PP2.coefs  = PP.coefs(n:(nt-1),:);


%-----------------------------------------------------------
% evaluates 2d spline at given S-values
function [XY] = myppval(PP, S);
[PPX, PPY] = one2two(PP);
X = ppval(PPX, S)';
Y = ppval(PPY, S)';
XY = [X,Y];


%-----------------------------------------------------------
% combines separate x,y splines into one 2d spline
function [PP] = two2one(PPX, PPY);

if ((PPX.pieces ~= PPY.pieces) | (PPX.order ~= PPY.order) | ...
      (PPX.dim ~= PPY.dim))
  fprintf('Error, PPX and PPY do not match.\n');
  return;
end

PP.breaks  = PPX.breaks;
PP.xcoefs  = PPX.coefs;
PP.ycoefs  = PPY.coefs;
PP.pieces  = PPX.pieces;
PP.order   = PPX.order;
PP.dim     = PPX.dim;

%-----------------------------------------------------------
% converts 2d spline into 2 1d splines
function [PPX, PPY] = one2two(PP);
PPX = mkpp(PP.breaks, PP.xcoefs);
PPY = mkpp(PP.breaks, PP.ycoefs);


%-----------------------------------------------------------
% Returns 10 1d quadrature points and weights
function [x, w] = quad1; 
x = zeros(10,1);
w = zeros(10,1);
x(0+1) = -0.9739065285171717200779640;
x(1+1) = -0.8650633666889845107320967;
x(2+1) = -0.6794095682990244062343274;
x(3+1) = -0.4333953941292471907992659;
x(4+1) = -0.1488743389816312108848260;
x(5+1) =  0.1488743389816312108848260 ;
x(6+1) =  0.4333953941292471907992659 ;
x(7+1) =  0.6794095682990244062343274 ;
x(8+1) =  0.8650633666889845107320967 ;
x(9+1) =  0.9739065285171717200779640 ;

w(0+1) =    0.0666713443086881375935688;
w(1+1) =   0.1494513491505805931457763 ;
w(2+1) =   0.2190863625159820439955349 ;
w(3+1) =   0.2692667193099963550912269 ;
w(4+1) =   0.2955242247147528701738930 ;
w(5+1) =  0.2955242247147528701738930  ;
w(6+1) =  0.2692667193099963550912269  ;
w(7+1) =  0.2190863625159820439955349  ;
w(8+1) =  0.1494513491505805931457763  ;
w(9+1) =  0.0666713443086881375935688  ;