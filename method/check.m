% check.m: Checks solutions found by solv.m for the real
%            business cycle model with indivisible labor
%            by making sure that the equilibrium conditions
%            hold in the steady state and after a technology
%            shock.  If the solution is correct, the program
%            will return a 6x2 vector of zeros.
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

% check steady state conditions

  sscond = ones(6,1);

  sscond(1) = gamma*css*hss - (1-theta)*yss;

  sscond(2) = eta - beta*(theta*(yss/kss)+1-delta);

  sscond(3) = yss - css - iss;

  sscond(4) = yss - a*(kss^theta)*(hss^(1-theta));

  sscond(5) = (eta-1+delta)*kss - iss;

  sscond(6) = a - a;

% construct period t values after shock

  et = sig;

  st = bigw*et;
  ft = bigu*st;

  kt = st(1);
  at = st(2);

  yt = ft(1);
  it = ft(2);
  ht = ft(3);
  ct = ft(4);

% construct period t+1 values after shock

  stp = bigpi*st;
  ftp = bigu*stp;

  ktp = stp(1);
  atp = stp(2);

  ytp = ftp(1);
  itp = ftp(2);
  htp = ftp(3);
  ctp = ftp(4);

% check equilibrium conditions

  eqcond = ones(6,1);

  eqcond(1) = ct + ht - yt;

  eqcond(2) = (eta/beta)*(ct-ctp) + (eta/beta-1+delta)*(ytp-ktp);

  eqcond(3) = (eta/beta-1+delta)*(yt-ct) + theta*(eta-1+delta)*(ct-it);

  eqcond(4) = yt - at - theta*kt - (1-theta)*ht;

  eqcond(5) = eta*ktp - (1-delta)*kt - (eta-1+delta)*it;

  eqcond(6) = atp - rho*at;

% report results

  [ sscond eqcond ]  