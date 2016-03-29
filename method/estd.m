% estd.m: Maximizes (minimizes negative) log likelihood function for
%           the real business cycle model with indivisible labor.
%
%         When maximizing the log likelihood function, the parameters
%           are transformed so that they satisfy theoretical 
%           restrictions.  The log likelihood function with
%           transformed parameters is contained in llfn.m.
%
%         The matrices D and V are constrained to be diagonal.
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

  global yt ct ht

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