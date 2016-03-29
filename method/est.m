% est.m: Maximizes (minimizes negative) log likelihood function for
%          the real business cycle model with indivisible labor.
%
%        When maximizing the log likelihood function, the parameters
%          are transformed so that they satisfy theoretical 
%          restrictions.  The log likelihood function with
%          transformed parameters is contained in llfn.m.
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

% load data and choose sample period

  global yt ct ht

  load ych.dat;

% full sample

  yt = ych(:,1);
  ct = ych(:,2);
  ht = ych(:,3);

% pre-1973 subsample

%  yt = ych(1:100,1);
%  ct = ych(1:100,2);
%  ht = ych(1:100,3);

% post-1973 subsample

%  yt = ych(101:218,1);
%  ct = ych(101:218,2);
%  ht = ych(101:218,3);

% pre-1980 subsample

%  yt = ych(1:128,1);
%  ct = ych(1:128,2);
%  ht = ych(1:128,3);

% post-1980 subsample

%  yt = ych(129:218,1);
%  ct = ych(129:218,2);
%  ht = ych(129:218,3);

% set starting values for full sample

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

% set starting values for pre-1973 sample

%  bettr = sqrt(0.99/(1-0.99));
%  gamtr = 0.0045;
%  thettr = sqrt(0.20/(1-0.20));
%  etatr = 0.0051;
%  delttr = sqrt(0.025/(1-0.025));
%  atr = 6;
%  rhotr = 0.9964;
%  sigtr = 0.0054;

%  dyytr = 1.1883;
%  dyctr = 0.3612;
%  dyhtr = -0.4075;
%  dcytr = 0.0875;
%  dcctr = 1.0024;
%  dchtr = -0.0738;
%  dhytr = 0.6527;
%  dhctr = 0.3161;
%  dhhtr = 0.3644;

%  vyytr = 0.0090;
%  vcytr = 0.0061;
%  vcctr = -0.0064;
%  vhytr = 0.0029;
%  vhctr = -0.0011;
%  vhhtr = -0.0031;

% set starting values for post-1973 sample

%  bettr = sqrt(0.99/(1-0.99));
%  gamtr = 0.0045;
%  thettr = sqrt(0.20/(1-0.20));
%  etatr = 0.0051;
%  delttr = sqrt(0.025/(1-0.025));
%  atr = 6;
%  rhotr = 0.9991;
%  sigtr = 0.0048;

%  dyytr = 1.2245;
%  dyctr = 0.7933;
%  dyhtr = -0.4947;
%  dcytr = 0.0285;
%  dcctr = 1.0665;
%  dchtr = -0.0560;
%  dhytr = 0.3008;
%  dhctr = 0.9692;
%  dhhtr = 0.3935;

%  vyytr = 0.0057;
%  vcytr = 0.0034;
%  vcctr = -0.0041;
%  vhytr = 0.0012;
%  vhctr = -0.0015;
%  vhhtr = 0.0000;

% set starting values for pre-1980 sample

%  bettr = sqrt(0.99/(1-0.99));
%  gamtr = 0.0045;
%  thettr = sqrt(0.20/(1-0.20));
%  etatr = 0.0051;
%  delttr = sqrt(0.025/(1-0.025));
%  atr = 6;
%  rhotr = 0.9959;
%  sigtr = 0.0052;

%  dyytr = 1.2635;
%  dyctr = 0.2386;
%  dyhtr = -0.4115;
%  dcytr = 0.0710;
%  dcctr = 1.0475;
%  dchtr = -0.1075;
%  dhytr = 0.6708;
%  dhctr = 0.2053;
%  dhhtr = 0.3779;

%  vyytr = -0.0087;
%  vcytr = -0.0056;
%  vcctr = 0.0065;
%  vhytr = -0.0030;
%  vhctr = 0.0006;
%  vhhtr = -0.0034;

% set starting values for post-1980 sample

%  bettr = sqrt(0.99/(1-0.99));
%  gamtr = 0.0045;
%  thettr = sqrt(0.20/(1-0.20));
%  etatr = 0.0051;
%  delttr = sqrt(0.025/(1-0.025));
%  atr = 6;
%  rhotr = 0.9985;
%  sigtr = 0.0041;

%  dyytr = 1.2801;
%  dyctr = 0.4877;
%  dyhtr = -0.4704;
%  dcytr = 0.0533;
%  dcctr = 1.0441;
%  dchtr = -0.0714;
%  dhytr = 0.3816;
%  dhctr = 0.6283;
%  dhhtr = 0.3914;

%  vyytr = -0.0062;
%  vcytr = -0.0026;
%  vcctr = -0.0042;
%  vhytr = -0.0016;
%  vhctr = -0.0007;
%  vhhtr = 0.0000;

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

%  thetstar = fminu('llfn2',bigtheto,options);

%  thetstar = fminu('llfn3',bigtheto,options);