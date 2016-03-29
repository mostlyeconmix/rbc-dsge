function llfn = llfn(bigthet);
% Uses the Kalman filter to evaluate the negative log likelihood
%   function for the real business cycle model with indivisible
%   labor.  Some parameters are transformed so that they satisty
%   theoretical restrictions.
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

% define variables and parameters

  global yt ct ht

  bigt = length(yt);

  bigthet = real(bigthet);

  bettr = sqrt(0.99/(1-0.99));
  gamtr = bigthet(1);
  thettr = bigthet(2);
  etatr = bigthet(3);
  delttr = sqrt(0.025/(1-0.025));
  atr = bigthet(4);
  rhotr = bigthet(5);
  sigtr = bigthet(6);

  dyytr = bigthet(7);
  dyctr = bigthet(8);
  dyhtr = bigthet(9);
  dcytr = bigthet(10);
  dcctr = bigthet(11);
  dchtr = bigthet(12);
  dhytr = bigthet(13);
  dhctr = bigthet(14);
  dhhtr = bigthet(15);

  vyytr = bigthet(16);
  vcytr = bigthet(17);
  vcctr = bigthet(18);
  vhytr = bigthet(19);
  vhctr = bigthet(20);
  vhhtr = bigthet(21);

% untransform parameters

  beta = bettr^2/(1+bettr^2);
  gamma = abs(gamtr);
  theta = thettr^2/(1+thettr^2);
  eta = 1 + abs(etatr);
  delta = delttr^2/(1+delttr^2);
  a = abs(atr);
  rho = rhotr;
  sig = abs(sigtr);

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

  bigdx = [ dyytr dyctr dyhtr ; ...
            dcytr dcctr dchtr ; ...
            dhytr dhctr dhhtr ];

  dxeig = eig(bigdx);

  dxviol = 0;

  if max(abs(dxeig)) > 1
    dxviol = 1;
  end

  bigv1x = sig^2;

  bigv2x1 = [ vyytr 0 0 ; ...
              vcytr vcctr 0 ; ...
              vhytr vhctr vhhtr ];

  bigv2x = bigv2x1*bigv2x1';

% form matrices FX, GX, and QX

  bigfx = [ bigax zeros(2,3) ; zeros(3,2) bigdx ];

  biggx = [ bigcx eye(3) ];

  bigqx = [ bigbx*bigv1x*bigbx' zeros(2,3) ; zeros(3,2) bigv2x ];

% put data in deviation form

  trend = 1:bigt;

  ythat = log(yt) - log(eta)*trend' - log(yss);

  cthat = log(ct) - log(eta)*trend' - log(css);

  hthat = log(ht) - log(hss);

  dthat = [ ythat cthat hthat ];

% evaluate negative log likelihood

  xt = zeros(5,1);

  bigsig1 = inv(eye(25)-kron(bigfx,bigfx))*bigqx(:);

  bigsigt = reshape(bigsig1,5,5);

  llfn = (3*bigt/2)*log(2*pi);

  for t = 1:bigt

    ut = dthat(t,:)' - biggx*xt;

    omegt = biggx*bigsigt*biggx';

    omeginvt = inv(omegt);

    llfn = llfn + (1/2)*(log(det(omegt))+ut'*omeginvt*ut);

    bigkt = bigfx*bigsigt*biggx'*omeginvt;

    xt = bigfx*xt + bigkt*ut;

    bigsigt = bigqx + bigfx*bigsigt*bigfx' ...
                - bigfx*bigsigt*biggx'*omeginvt*biggx*bigsigt*bigfx';

  end

% penalize constraint violations

  if dxviol

    llfn = llfn + 1e8;

  end

  if abs(imag(llfn)) > 0

    llfn = real(llfn) + 1e8;

  end

  if abs(rhotr) > 1

    llfn = llfn + 1e8;

  end
