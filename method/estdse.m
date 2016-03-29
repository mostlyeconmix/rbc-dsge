% estdse.m: Maximizes (minimizes negative) log likelihood function for
%            the real business cycle model with indivisible labor
%            and computes standard errors for estimated parameters.
%
%           When maximizing the log likelihood function, the parameters
%            are transformed so that they satisfy theoretical 
%            restrictions.  The log likelihood function with
%            transformed parameters is contained in llfn.m.
%
%           When calculating standard errors, the original
%            parameters are used.  The log likelihood function without
%            transformed parameters is contained in llfnse.m.
%
%           The untransformed parameter estimates are reported as
%            the vector tstar.  The standard errors are reported as
%            the vector sevec.
%
%           The matrices D and V are constrained to be diagonal.
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

% load data

  global yt ct ht scalinv

  load ych.dat;

  yt = ych(:,1);
  ct = ych(:,2);
  ht = ych(:,3);

% set starting values

  bettr = sqrt(0.99/(1-0.99));
  gamtr = 0.0045;
  thettr = sqrt(0.20/(1-0.20));
  etatr = 0.0051;
  delttr = sqrt(0.025/(1-0.025));
  atr = 6;
  rhotr = 0.9914;
  sigtr = 0.0073;

  dyytr = -0.9845;
  dyctr = 0;
  dyhtr = 0;
  dcytr = 0;
  dcctr = 0.9836;
  dchtr = 0;
  dhytr = 0;
  dhctr = 0;
  dhhtr = 0.9935;

  vyytr = 0.0001;
  vcytr = 0;
  vcctr = 0.0061;
  vhytr = 0;
  vhctr = 0;
  vhhtr = 0.0078;

% maximize likelihood

  bigtheto = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dcctr dhhtr ...
               vyytr vcctr vhhtr ]';

  options(1) = 1;
  options(14) = 10000;

  thetstar = fminu('llfnd',bigtheto,options);

% find standard errors

  thetstar = real(thetstar);

  beta = bettr^2/(1+bettr^2);
  gamma = abs(thetstar(1));
  theta = thetstar(2)^2/(1+thetstar(2)^2);
  eta = 1 + abs(thetstar(3));
  delta = delttr^2/(1+delttr^2);
  a = abs(thetstar(4));
  rho = thetstar(5);
  sig = abs(thetstar(6));

  dyy = thetstar(7);
  dyc = 0;
  dyh = 0;
  dcy = 0;
  dcc = thetstar(8);
  dch = 0;
  dhy = 0;
  dhc = 0;
  dhh = thetstar(9);

  cholv = [ thetstar(10) 0 0 ; ...
            0 thetstar(11) 0 ; ...
            0 0 thetstar(12) ];

  bigv = cholv*cholv';

  vyy = sqrt(bigv(1,1));
  vcc = sqrt(bigv(2,2));
  vhh = sqrt(bigv(3,3));
  vyc = bigv(1,2);
  vyh = bigv(1,3);
  vch = bigv(2,3);

  tstar = [ gamma theta eta a rho sig ...
            dyy dcc dhh ...
            vyy vcc vhh ]';

  scalvec = ones(12,1);

  scalinv = inv(diag(scalvec));

  tstars = diag(scalvec)*tstar;

  fstar = llfndse(tstars);

  eee = 1e-6;

  epsmat = eee*eye(12);

  hessvec = zeros(12,1);

  for i = 1:12

    hessvec(i) = llfndse(tstars+epsmat(:,i));

  end

  hessmat = zeros(12,12);

  for i = 1:12

    for j = 1:12

      hessmat(i,j) = (llfndse(tstars+epsmat(:,i)+epsmat(:,j)) ...
                        -hessvec(i)-hessvec(j)+fstar)/eee^2;

    end

  end

  bighx = scalinv*inv(hessmat)*scalinv';

  sevec = sqrt(diag(bighx));