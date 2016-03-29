% estseqd.m: Estimates the parameters of the real business cycle model
%              with indivisible labor for a sequence of sample periods
%              of increasing length.  Starts with the sample period
%              1948:1-1984:4 and ends with the sample period
%              1948:1-2002:2.  Thus, estimates the parameters for a
%              total of 71 samples.
%
%            The matrices D and V are constrained to be diagonal.
%
%            Returns:
%              thetmatd = columns contain transformed estimates for each
%                                                                  sample
%              tstrmatd = columns contain untransformed estimates for
%                                                             each sample
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

% load data and set up output matrix

  global yt ct ht

  load ych.dat;

  bigt = length(ych(:,1));

  thetmatd = zeros(12,71);

  tstrmatd = zeros(21,71);

% set starting values

  gamtr = 0.0045;
  thettr = sqrt(0.20/(1-0.20));
  etatr = 0.0051;
  atr = 6;
  rhotr = 0.9914;
  sigtr = 0.0073;

  dyytr = -0.9845;
  dcctr = 0.9836;
  dhhtr = 0.9935;

  vyytr = 0.0001;
  vcctr = 0.0061;
  vhhtr = 0.0078;

  bigthet1 = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dcctr dhhtr vyytr vcctr vhhtr ]';

  gamtr = 0.0039;
  thettr = 0.9693;
  etatr = 0.0049;
  atr = 1.1789;
  rhotr = 0.9987;
  sigtr = 0.0093;

  dyytr = -0.9831;
  dcctr = 0.9997;
  dhhtr = 0.9892;

  vyytr = 0.0001;
  vcctr = 0.0063;
  vhhtr = 0.0076;

  bigthet2 = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dcctr dhhtr vyytr vcctr vhhtr ]';

  gamtr = 0.0041;
  thettr = 0.9157;
  etatr = 0.0051;
  atr = 1.4340;
  rhotr = 0.9981;
  sigtr = 0.0099;

  dyytr = -0.9934;
  dcctr = 0.9996;
  dhhtr = 0.9576;

  vyytr = 0.0001;
  vcctr = 0.0068;
  vhhtr = 0.0083;

  bigthet3 = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dcctr dhhtr vyytr vcctr vhhtr ]';

  gamtr = 0.0041;
  thettr = 0.9092;
  etatr = 0.0050;
  atr = 1.4646;
  rhotr = 0.9981;
  sigtr = 0.0100;

  dyytr = -0.9999;
  dcctr = 0.9995;
  dhhtr = 0.9560;

  vyytr = 0.0001;
  vcctr = 0.0068;
  vhhtr = 0.0083;

  bigthet4 = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dcctr dhhtr vyytr vcctr vhhtr ]';

% estimate parameters for each sample period

  for t = bigt:-1:148

    t

    yt = ych(1:t,1);
    ct = ych(1:t,2);
    ht = ych(1:t,3);

    options(1) = 1;
    options(14) = 100000;

    if t == bigt

      thetstar = fminu('llfnd',bigthet1,options);

    elseif t == 217

      thetstar = fminu('llfnd',bigthet1,options);

    elseif t == 214

      thetstar = fminu('llfnd',bigthet1,options);

    elseif t == 192

      thetstar = fminu('llfnd',bigthet2,options);

    elseif t == 186

      thetstar = fminu('llfnd',bigthet2,options);

    elseif t == 151

      thetstar = fminu('llfnd',bigthet3,options);

    elseif t == 148

      thetstar = fminu('llfnd',bigthet4,options);

    else

      thetstar = fminu('llfnd',thetstar,options);

    end

    thetstar = real(thetstar);

    gamma = abs(thetstar(1));
    theta = thetstar(2)^2/(1+thetstar(2)^2);
    eta = 1 + abs(thetstar(3));
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
              dyy dyc dyh dcy dcc dch dhy dhc dhh ...
              vyy vcc vhh vyc vyh vch ]';

    thetmatd(:,t-147) = thetstar;

    tstrmatd(:,t-147) = tstar;

  end

% save results

  save thetmatd thetmatd

  save tstrmatd tstrmatd