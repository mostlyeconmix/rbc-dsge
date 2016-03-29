% stabtest.m: Performs parameter stability tests for the real
%               business cycle model with indivisible labor.
%               Breakpoint occurs at 1973:1.
%
%             Returns:
%
%               tstar1 = estimated parameters, 1948:1-1972:4
%               sevec1 = standard errors, 1948:1-1972:4
%               tstar2 = estimated parameters, 1973:1-2002:2
%               sevec2 = standard errors, 1973:1-2002:2
%               lrstat = LR statistic for stability of all parameters
%               wstat = Wald statistic for stability of all parameters
%               wstatq = Wald statistic for stability of structural
%                                                         parameters
%               wstatx = Wald statistic for stability of structural
%                          parameters, excluding average productivity growth
%               wstatn = Wald statistic for stability of other
%                                                         parameters
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

% define pre-1973 subsample

  yt = yt(1:100);
  ct = ct(1:100);
  ht = ht(1:100);

% set starting values for pre-1973 subsample

  bettr = sqrt(0.99/(1-0.99));
  gamtr = 0.0045;
  thettr = sqrt(0.20/(1-0.20));
  etatr = 0.0051;
  delttr = sqrt(0.025/(1-0.025));
  atr = 6;
  rhotr = 0.9964;
  sigtr = 0.0054;

  dyytr = 1.1883;
  dyctr = 0.3612;
  dyhtr = -0.4075;
  dcytr = 0.0875;
  dcctr = 1.0024;
  dchtr = -0.0738;
  dhytr = 0.6527;
  dhctr = 0.3161;
  dhhtr = 0.3644;

  vyytr = 0.0090;
  vcytr = 0.0061;
  vcctr = -0.0064;
  vhytr = 0.0029;
  vhctr = -0.0011;
  vhhtr = -0.0031;

% maximize likelihood for pre-1973 subsample

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

% find standard errors for pre-1973 subsample

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

  tstar1 = [ gamma theta eta a rho sig ...
             dyy dyc dyh dcy dcc dch dhy dhc dhh ...
             vyy vcc vhh vyc vyh vch ]';

  scalvec = [ ones(19,1) ; 10 ; 10 ];

  scalinv = inv(diag(scalvec));

  tstar1s = diag(scalvec)*tstar1;

  fstar1 = llfnse(tstar1s);

  eee = 1e-6;

  epsmat = eee*eye(21);

  hessvec = zeros(21,1);

  for i = 1:21

    hessvec(i) = llfnse(tstar1s+epsmat(:,i));

  end

  hessmat = zeros(21,21);

  for i = 1:21

    for j = 1:21

      hessmat(i,j) = (llfnse(tstar1s+epsmat(:,i)+epsmat(:,j)) ...
                        -hessvec(i)-hessvec(j)+fstar1)/eee^2;

    end

  end

  bighx1 = scalinv*inv(hessmat)*scalinv';

  sevec1 = sqrt(diag(bighx1));

% reload data

  load ych.dat;

  yt = ych(:,1);
  ct = ych(:,2);
  ht = ych(:,3);

% define post-1973 subsample

  yt = yt(101:218);
  ct = ct(101:218);
  ht = ht(101:218);

% set starting values for post-1973 subsample

  bettr = sqrt(0.99/(1-0.99));
  gamtr = 0.0045;
  thettr = sqrt(0.20/(1-0.20));
  etatr = 0.0051;
  delttr = sqrt(0.025/(1-0.025));
  atr = 6;
  rhotr = 0.9991;
  sigtr = 0.0048;

  dyytr = 1.2245;
  dyctr = 0.7933;
  dyhtr = -0.4947;
  dcytr = 0.0285;
  dcctr = 1.0665;
  dchtr = -0.0560;
  dhytr = 0.3008;
  dhctr = 0.9692;
  dhhtr = 0.3935;

  vyytr = 0.0057;
  vcytr = 0.0034;
  vcctr = -0.0041;
  vhytr = 0.0012;
  vhctr = -0.0015;
  vhhtr = 0.0000;

% maximize likelihood for post-1973 subsample

  bigtheto = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dyctr dyhtr ...
               dcytr dcctr dchtr ...
               dhytr dhctr dhhtr ...
               vyytr  ...
               vcytr vcctr ...
               vhytr vhctr vhhtr ]';

  options(1) = 1;
  options(14) = 10000;

  thetstar = fminu('llfn2',bigtheto,options);

% find standard errors for post-1973 subsample

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

  tstar2 = [ gamma theta eta a rho sig ...
             dyy dyc dyh dcy dcc dch dhy dhc dhh ...
             vyy vcc vhh vyc vyh vch ]';

  scalvec = ones(21,1);

  scalinv = inv(diag(scalvec));

  tstar2s = diag(scalvec)*tstar2;

  fstar2 = llfnse2(tstar2);

  eee = 1e-6;

  epsmat = eee*eye(21);

  hessvec = zeros(21,1);

  for i = 1:21

    hessvec(i) = llfnse2(tstar2s+epsmat(:,i));

  end

  hessmat = zeros(21,21);

  for i = 1:21

    for j = 1:21

      hessmat(i,j) = (llfnse2(tstar2s+epsmat(:,i)+epsmat(:,j)) ...
                        -hessvec(i)-hessvec(j)+fstar2)/eee^2;

    end

  end

  bighx2 = scalinv*inv(hessmat)*scalinv';

  sevec2 = sqrt(diag(bighx2));

% form test statistics

  lrstat = -2*(fstar1+fstar2+2323.5501);

  wstat = (tstar1-tstar2)'*inv(bighx1+bighx2)*(tstar1-tstar2);

  bighxq1 = bighx1(1:6,1:6);
  bighxq2 = bighx2(1:6,1:6);

  wstatq = (tstar1(1:6)-tstar2(1:6))' ...
             *inv(bighxq1+bighxq2)*(tstar1(1:6)-tstar2(1:6));
        
  bighxx1 = [ bighx1(1:2,1:2) bighx1(1:2,4:6) ; ...
              bighx1(4:6,1:2) bighx1(4:6,4:6) ];
      
  bighxx2 = [ bighx2(1:2,1:2) bighx2(1:2,4:6) ; ...
              bighx2(4:6,1:2) bighx2(4:6,4:6) ];    
  
  tstarx1 = [ tstar1(1:2) ; tstar1(4:6) ];
  
  tstarx2 = [ tstar2(1:2) ; tstar2(4:6) ];

  wstatx = (tstarx1-tstarx2)' ...
             *inv(bighxx1+bighxx2)*(tstarx1-tstarx2);

  bighxn1 = bighx1(7:21,7:21);
  bighxn2 = bighx2(7:21,7:21);

  wstatn = (tstar1(7:21)-tstar2(7:21))' ...
             *inv(bighxn1+bighxn2)*(tstar1(7:21)-tstar2(7:21));