% estse.m: Maximizes (minimizes negative) log likelihood function for
%           the real business cycle model with indivisible labor
%           and computes standard errors for estimated parameters.
%
%         When maximizing the log likelihood function, the parameters
%           are transformed so that they satisfy theoretical 
%           restrictions.  The log likelihood function with
%           transformed parameters is contained in llfn.m.
%
%         When calculating standard errors, the original
%            parameters are used.  The log likelihood function without
%            transformed parameters is contained in llfnse.m.
%
%         The untransformed parameter estimates are reported as
%            the vector tstar.  The standard errors are reported as
%            the vector sevec.
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
  rhotr = 0.9975;
  sigtr = 0.0055;

  dyytr = 1.4187;
  dyctr = 0.2251;
  dyhtr = -0.4441;
  dcytr = 0.0935;
  dcctr = 1.0236;
  dchtr = -0.0908;
  dhytr = 0.7775;
  dhctr = 0.3706;
  dhhtr = 0.2398;

  vyytr = 0.0072;
  vcytr = 0.0040;
  vcctr = 0.0057;
  vhytr = 0.0015;
  vhctr = 0.0010;
  vhhtr = 0.0000;

% maximize likelihood

  bigtheto = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dyctr dyhtr ...
               dcytr dcctr dchtr ...
               dhytr dhctr dhhtr ...
               vyytr  ...
               vcytr vcctr ...
               vhytr vhctr vhhtr ]';

  options(1) = 1;
  options(14) = 10000;

  thetstar = fminu('llfn',bigtheto,options);

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
  dyc = thetstar(8);
  dyh = thetstar(9);
  dcy = thetstar(10);
  dcc = thetstar(11);
  dch = thetstar(12);
  dhy = thetstar(13);
  dhc = thetstar(14);
  dhh = thetstar(15);

  cholv = [ thetstar(16) 0 0 ; ...
            thetstar(17:18)' 0 ; ...
            thetstar(19:21)' ];

  bigv = cholv*cholv';

  vyy = sqrt(bigv(1,1));
  vcc = sqrt(bigv(2,2));
  vhh = sqrt(bigv(3,3));
  vyc = bigv(1,2);
  vyh = bigv(1,3);
  vch = bigv(2,3);

  tstar = [ gamma theta eta a rho sig ...
            dyy dyc dyh dcy dcc dch dhy dhc dhh ...
            vyy vcc vhh vyc vyh vch ]';

  scalvec = ones(21,1);

  scalinv = inv(diag(scalvec));

  tstars = diag(scalvec)*tstar;

  fstar = llfnse(tstars);

  eee = 1e-6;

  epsmat = eee*eye(21);

  hessvec = zeros(21,1);

  for i = 1:21

    hessvec(i) = llfnse(tstars+epsmat(:,i));

  end

  hessmat = zeros(21,21);

  for i = 1:21

    for j = 1:21

      hessmat(i,j) = (llfnse(tstars+epsmat(:,i)+epsmat(:,j)) ...
                        -hessvec(i)-hessvec(j)+fstar)/eee^2;

    end

  end

  bighx = scalinv*inv(hessmat)*scalinv';

  sevec = sqrt(diag(bighx));