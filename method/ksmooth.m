% ksmooth.m: Uses the Kalman smoother, together with the untransformed estimates
%              found be estse.m, to produce a series of smoothed estimates
%              of the hybrid model's shocks.
%
%            The results are contained in the 5xT matrix xbigtvec, which
%              keeps track of the five state variables in the hybrid model's
%              state equation.
%
%            Also calculates corey, corec, coreh, the correlation
%              coefficients between the innovation to technology and the 
%              innovations to each of the three residuals.
%
% THIS PROGRAM WAS WRITTEN FOR MATLAB BY
%
%   PETER N. IRELAND
%   BOSTON COLLEGE
%   DEPARTMENT OF ECONOMICS
%   140 COMMONWEALTH AVENUE
%   CHESTNUT HILL, MA 02467
%   irelandp@bc.edu
%
%  FINANCIAL SUPPORT FROM THE NATIONAL SCIENCE FOUNDATION UNDER
%    GRANT NOS. SES-9985763 AND SES-0213461 IS GRATEFULLY ACKNOWLEDGED.
%
%  COPYRIGHT (c) 2003 BY PETER N. IRELAND.  REDISTRIBUTION IS
%    PERMITTED FOR EDUCATIONAL AND RESEARCH PURPOSES, SO LONG AS
%    NO CHANGES ARE MADE.  ALL COPIES MUST BE PROVIDED FREE OF
%    CHARGE AND MUST INCLUDE THIS COPYRIGHT NOTICE.

% compute steady-state values

  kappa = eta/beta - 1 + delta;

  lambda = eta - 1 + delta;

  hss = ((1-theta)/gamma)/(1-theta*lambda/kappa);

  yss = (a^(1/(1-theta)))*((theta/kappa)^(theta/(1-theta)))*hss;

  kss = (theta/kappa)*yss;

  iss = (theta*lambda/kappa)*yss;

  css = yss - iss;

% compute K coefficients

  bigk11 = (eta-beta*(1-theta)*(1-delta))/(beta*eta*theta);

  bigk12 = (beta*eta*theta^2-eta+beta*(1-theta^2)*(1-delta)) ...
             /(beta*eta*theta^2);

  bigk22 = eta*theta/(eta-beta*(1-theta)*(1-delta));

% compute L coefficients

  bigl1 = (eta-beta*(1-delta))/(beta*eta*theta^2);

  bigl2 = rho*(eta-beta*(1-delta))/(eta-beta*(1-theta)*(1-delta));

% form S matrices

  bigs1 = (bigk22-bigk11)/bigk12;

  bigs2 = ((bigk22-bigk11)*bigl1-bigk12*bigl2)/(bigk12*(bigk11-rho));

  bigs3 = bigk22;

  bigs4 = bigk12*bigl2/(bigk22-bigk11) ...
            + (bigk22-rho)*((bigk22-bigk11)*bigl1-bigk12*bigl2) ...
            / ((bigk22-bigk11)*(bigk11-rho));

  bigs5 = [ 1 - ((1-theta)/theta)*bigs1 ; ...
            bigs1 + (kappa/(theta^2*lambda))*(theta-bigs1) ; ...
            1 - (1/theta)*bigs1 ];

  bigs6 = [ 1/theta - ((1-theta)/theta)*bigs2 ; ...
            bigs2 + (kappa/(theta^2*lambda))*(1-bigs2) ; ...
            1/theta - (1/theta)*bigs2 ];

% form matrices PI, W, and U

  bigpi = [ bigs3 bigs4 ; 0 rho ];

  bigw = [ 0 ; 1 ];

  bigu = [ bigs5 bigs6 ; bigs1 bigs2 ];

% form matrices AX, BX, CX, DX, V1X, and V2x

  bigax = bigpi;

  bigbx = bigw;

  bigcx = [ bigu(1,:) ; bigu(4,:) ; bigu(3,:) ];

  bigdx = [ dyy dyc dyh ; ...
            dcy dcc dch ; ...
            dhy dhc dhh ];

  dxeig = eig(bigdx);

  dxviol = 0;

  if max(abs(dxeig)) > 1
    dxviol = 1;
  end

  bigv1x = sig^2;
  
  bigv2x = [ vyy^2 vyc vyh ; ...
             vyc vcc^2 vch ; ...
             vyh vch vhh^2 ];

% form matrices FX, GX, and QX

  bigfx = [ bigax zeros(2,3) ; zeros(3,2) bigdx ];

  biggx = [ bigcx eye(3) ];

  bigqx = [ bigbx*bigv1x*bigbx' zeros(2,3) ; zeros(3,2) bigv2x ];

% put data in deviation form

  bigt = length(yt);

  trend = 1:bigt;

  ythat = log(yt) - log(eta)*trend' - log(yss);

  cthat = log(ct) - log(eta)*trend' - log(css);

  hthat = log(ht) - log(hss);

  dthat = [ ythat cthat hthat ];

% run through the Kalman Filter
  
  xttvec = zeros(5,bigt);
  
  xtlvec = zeros(5,bigt);
  
  sigttvec = zeros(25,bigt);
 
  sigtlvec = zeros(25,bigt);

  xt = zeros(5,1);

  bigsigt = inv(eye(25)-kron(bigfx,bigfx))*bigqx(:);
  
  bigsigt = reshape(bigsigt,5,5);
  
  for t = 1:bigt
      
    xtlvec(:,t) = xt;
    
    sigtlvec(:,t) = bigsigt(:);
      
    omegt = biggx*bigsigt*biggx';
    
    omeginvt = inv(omegt);
    
    ut = dthat(t,:)' - biggx*xt;
    
    xtt = xt + bigsigt*biggx'*omeginvt*ut;
    
    bigsigtt = bigsigt - bigsigt*biggx'*omeginvt*biggx*bigsigt;
    
    xttvec(:,t) = xtt;
    
    sigttvec(:,t) = bigsigtt(:);
    
    xt = bigfx*xtt;
    
    bigsigt = bigqx + bigfx*bigsigtt*bigfx';
    
  end  

% run through Kalman smoother

  xbigtvec = zeros(5,bigt);
  
  xbigtvec(:,bigt) = xttvec(:,bigt);
  
  for j = 1:bigt-1
      
    bigsigtt = sigttvec(:,bigt-j);
    
    bigsigtt = reshape(bigsigtt,5,5);
    
    bigsigtp = sigtlvec(:,bigt-j+1);
    
    bigsigtp = reshape(bigsigtp,5,5);
    
    bigjt = bigsigtt*bigfx'*inv(bigsigtp);
    
    xtt = xttvec(:,bigt-j);
    
    xtpbigt = xbigtvec(:,bigt-j+1);
    
    xtpt = xtlvec(:,bigt-j+1);
    
    xbigtvec(:,bigt-j)= xtt + bigjt*(xtpbigt-xtpt);
    
  end  
  
% calculate correlations between innovations

  etavec = xbigtvec(:,2:bigt)' - xbigtvec(:,1:bigt-1)'*bigfx';
  
  epsvec = etavec(:,2);
  
  xiyvec = etavec(:,3);
  
  xicvec = etavec(:,4);
  
  xihvec = etavec(:,5);
  
  corey = corrcoef(epsvec,xiyvec);
  
  corec = corrcoef(epsvec,xicvec);
  
  coreh = corrcoef(epsvec,xihvec);