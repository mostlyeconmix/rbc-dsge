% estseq.m: Estimates the parameters of the real business cycle model
%             with indivisible labor for a sequence of sample periods
%             of increasing length.  Starts with the sample period
%             1948:1-1984:4 and ends with the sample period
%             1948:1-2002:2.  Thus, estimates the parameters for a
%             total of 71 samples.
%
%           Returns:
%             thetmat = columns contain transformed estimates for each
%                                                                 sample
%             tstarmat = columns contain untransformed estimates for
%                                                            each sample
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

  thetmat = zeros(21,71);

  tstarmat = zeros(21,71);

% set starting values

  gamtr = 0.0045;
  thettr = sqrt(0.20/(1-0.20));
  etatr = 0.0051;
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

  bigthet1 = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dyctr dyhtr ...
               dcytr dcctr dchtr ...
               dhytr dhctr dhhtr ...
               vyytr  ...
               vcytr vcctr ...
               vhytr vhctr vhhtr ]';

  gamtr = 0.0045;
  thettr = 0.5441;
  etatr = 0.0050;
  atr = 5.2174;
  rhotr = 0.9986;
  sigtr = 0.0057;

  dyytr = 1.3897;
  dyctr = 0.3617;
  dyhtr = -0.5032;
  dcytr = 0.1471;
  dcctr = 0.9692;
  dchtr = -0.1138;
  dhytr = 0.7510;
  dhctr = 0.4285;
  dhhtr = 0.2050;

  vyytr = 0.0071;
  vcytr = 0.0044;
  vcctr = 0.0056;
  vhytr = 0.0014;
  vhctr = 0.0012;
  vhhtr = 0.0000;

  bigthet2 = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dyctr dyhtr ...
               dcytr dcctr dchtr ...
               dhytr dhctr dhhtr ...
               vyytr  ...
               vcytr vcctr ...
               vhytr vhctr vhhtr ]';

  gamtr = 0.0046;
  thettr = 0.5367;
  etatr = 0.0048;
  atr = 5.3634;
  rhotr = 0.9987;
  sigtr = 0.0058;

  dyytr = 1.3672;
  dyctr = 0.3997;
  dyhtr = -0.5324;
  dcytr = 0.1404;
  dcctr = 0.9805;
  dchtr = -0.1230;
  dhytr = 0.7083;
  dhctr = 0.4577;
  dhhtr = 0.2097;

  vyytr = 0.0070;
  vcytr = 0.0046;
  vcctr = 0.0056;
  vhytr = 0.0012;
  vhctr = 0.0014;
  vhhtr = 0.0000;

  bigthet3 = [ gamtr thettr etatr atr rhotr sigtr ...
               dyytr dyctr dyhtr ...
               dcytr dcctr dchtr ...
               dhytr dhctr dhhtr ...
               vyytr  ...
               vcytr vcctr ...
               vhytr vhctr vhhtr ]';

% estimate parameters for each sample period

  for t = bigt:-1:148

    t

    yt = ych(1:t,1);
    ct = ych(1:t,2);
    ht = ych(1:t,3);

    options(1) = 1;
    options(14) = 100000;

    if t == bigt

      thetstar = fminu('llfn',bigthet1,options);

    elseif t == 217

      thetstar = fminu('llfn',bigthet1,options);

    elseif t == 210

      thetstar = fminu('llfn',bigthet1,options);

    elseif t == 207

      thetstar = fminu('llfn',bigthet2,options);

    elseif t == 206

      thetstar = fminu('llfn',bigthet2,options);

    elseif t == 205

      thetstar = fminu('llfn',bigthet2,options);

    elseif t == 201

      thetstar = fminu('llfn',bigthet2,options);

    elseif t == 199

      thetstar = fminu('llfn',bigthet2,options);

    elseif t == 195

      thetstar = fminu('llfn',bigthet2,options);

    elseif t == 192

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 191

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 188

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 186

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 183

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 179

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 172

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 158

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 157

      thetstar = fminu('llfn',bigthet3,options);

    elseif t == 152

      thetstar = fminu('llfn',bigthet3,options);

    else

      thetstar = fminu('llfn',thetstar,options);

    end

    thetstar = real(thetstar);

    gamma = abs(thetstar(1));
    theta = thetstar(2)^2/(1+thetstar(2)^2);
    eta = 1 + abs(thetstar(3));
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

    thetmat(:,t-147) = thetstar;

    tstarmat(:,t-147) = tstar;

  end

% save results

  save thetmat thetmat

  save tstarmat tstarmat