function vardecfn = vardecfn(bigthet);
% For given parameter estimates, computes variance decompositions
%   for the real business cycle model with indivisible labor.
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

% define parameters

  beta = 0.99;
  gamma = bigthet(1);
  theta = bigthet(2);
  eta = bigthet(3);
  delta = 0.025;
  a = bigthet(4);
  rho = bigthet(5);
  sig = bigthet(6);

  dyy = bigthet(7);
  dyc = bigthet(8);
  dyh = bigthet(9);
  dcy = bigthet(10);
  dcc = bigthet(11);
  dch = bigthet(12);
  dhy = bigthet(13);
  dhc = bigthet(14);
  dhh = bigthet(15);

  vyy = bigthet(16);
  vcc = bigthet(17);
  vhh = bigthet(18);
  vyc = bigthet(19);
  vyh = bigthet(20);
  vch = bigthet(21);

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

% calculate k-step ahead forecast error variances

  vardecfn = zeros(56,1);

  bigf2x = eye(5);

  bigg2x = zeros(4,5);
  bigg2x(1:2,:) = biggx(1:2,:);
  bigg2x(3,:) = (yss/iss)*biggx(1,:) - (css/iss)*biggx(2,:);
  bigg2x(4,:) = biggx(3,:);

  bigsigxk = zeros(5,5);

  for k = 1:40

    bigsigxk = bigsigxk + bigf2x*bigqx*bigf2x';

    bigsigyk = bigg2x*bigsigxk*bigg2x';

    if k == 1
      vardecfn(1) = bigsigyk(1,1);
      vardecfn(8) = bigsigyk(2,2);
      vardecfn(15) = bigsigyk(3,3);
      vardecfn(22) = bigsigyk(4,4);
    end

    if k == 4
      vardecfn(2) = bigsigyk(1,1);
      vardecfn(9) = bigsigyk(2,2);
      vardecfn(16) = bigsigyk(3,3);
      vardecfn(23) = bigsigyk(4,4);
    end

    if k == 8
      vardecfn(3) = bigsigyk(1,1);
      vardecfn(10) = bigsigyk(2,2);
      vardecfn(17) = bigsigyk(3,3);
      vardecfn(24) = bigsigyk(4,4);
    end

    if k == 12
      vardecfn(4) = bigsigyk(1,1);
      vardecfn(11) = bigsigyk(2,2);
      vardecfn(18) = bigsigyk(3,3);
      vardecfn(25) = bigsigyk(4,4);
    end

    if k == 20
      vardecfn(5) = bigsigyk(1,1);
      vardecfn(12) = bigsigyk(2,2);
      vardecfn(19) = bigsigyk(3,3);
      vardecfn(26) = bigsigyk(4,4);
    end

    if k == 40
      vardecfn(6) = bigsigyk(1,1);
      vardecfn(13) = bigsigyk(2,2);
      vardecfn(20) = bigsigyk(3,3);
      vardecfn(27) = bigsigyk(4,4);
    end

    bigf2x = bigf2x*bigfx;

  end    

  bigsigx1 = inv(eye(25)-kron(bigfx,bigfx))*bigqx(:);

  bigsigx = reshape(bigsigx1,5,5);

  bigsigy = bigg2x*bigsigx*bigg2x';

  vardecfn(7) = bigsigy(1,1);
  vardecfn(14) = bigsigy(2,2);
  vardecfn(21) = bigsigy(3,3);
  vardecfn(28) = bigsigy(4,4);

% calculate fractions due to technology

  bigqxt = [ bigbx*bigv1x*bigbx' zeros(2,3) ; zeros(3,5) ];

  bigf2x = eye(5);

  bigsigxk = zeros(5,5);

  for k = 1:40

    bigsigxk = bigsigxk + bigf2x*bigqxt*bigf2x';

    bigsigyk = bigg2x*bigsigxk*bigg2x';

    if k == 1
      vardecfn(29) = bigsigyk(1,1)/vardecfn(1);
      vardecfn(36) = bigsigyk(2,2)/vardecfn(8);
      vardecfn(43) = bigsigyk(3,3)/vardecfn(15);
      vardecfn(50) = bigsigyk(4,4)/vardecfn(22);
    end

    if k == 4
      vardecfn(30) = bigsigyk(1,1)/vardecfn(2);
      vardecfn(37) = bigsigyk(2,2)/vardecfn(9);
      vardecfn(44) = bigsigyk(3,3)/vardecfn(16);
      vardecfn(51) = bigsigyk(4,4)/vardecfn(23);
    end

    if k == 8
      vardecfn(31) = bigsigyk(1,1)/vardecfn(3);
      vardecfn(38) = bigsigyk(2,2)/vardecfn(10);
      vardecfn(45) = bigsigyk(3,3)/vardecfn(17);
      vardecfn(52) = bigsigyk(4,4)/vardecfn(24);
    end

    if k == 12
      vardecfn(32) = bigsigyk(1,1)/vardecfn(4);
      vardecfn(39) = bigsigyk(2,2)/vardecfn(11);
      vardecfn(46) = bigsigyk(3,3)/vardecfn(18);
      vardecfn(53) = bigsigyk(4,4)/vardecfn(25);
    end

    if k == 20
      vardecfn(33) = bigsigyk(1,1)/vardecfn(5);
      vardecfn(40) = bigsigyk(2,2)/vardecfn(12);
      vardecfn(47) = bigsigyk(3,3)/vardecfn(19);
      vardecfn(54) = bigsigyk(4,4)/vardecfn(26);
    end

    if k == 40
      vardecfn(34) = bigsigyk(1,1)/vardecfn(6);
      vardecfn(41) = bigsigyk(2,2)/vardecfn(13);
      vardecfn(48) = bigsigyk(3,3)/vardecfn(20);
      vardecfn(55) = bigsigyk(4,4)/vardecfn(27);
    end

    bigf2x = bigf2x*bigfx;

  end    

  bigsigx1 = inv(eye(25)-kron(bigfx,bigfx))*bigqxt(:);

  bigsigx = reshape(bigsigx1,5,5);

  bigsigy = bigg2x*bigsigx*bigg2x';

  vardecfn(35) = bigsigy(1,1)/vardecfn(7);
  vardecfn(42) = bigsigy(2,2)/vardecfn(14);
  vardecfn(49) = bigsigy(3,3)/vardecfn(21);
  vardecfn(56) = bigsigy(4,4)/vardecfn(28);

% convert to percentages

  vardecfn = 100*vardecfn;