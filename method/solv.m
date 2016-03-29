% solv.m: Solves the real business cycle model with indivisible labor.
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

% set parameter values

  beta = 0.99;
  gamma = 0.0045;
  theta = 0.20;
  eta = 1.0042;
  delta = 0.025;
  a = 6;
  rho = 0.95;
  sig = 0.007;

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