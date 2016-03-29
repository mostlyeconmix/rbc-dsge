% imp.m: Computes impulse responses for the real business
%          cycle model with indivisible labor using the
%          solutions found by solv.m.
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

% technology shock

  st = zeros(50,2);
  ft = zeros(50,4);

  et = sig;

  st(2,:) = (bigw*et)';
  ft(2,:) = (bigu*st(2,:)')';

  for t = 3:50;

    st(t,:) = (bigpi*st(t-1,:)')';
    ft(t,:) = (bigu*st(t,:)')';

  end

% create output vectors

  st = 100*st;
  ft = 100*ft;

  kt = st(:,1);
  at = st(:,2);

  yt = ft(:,1);
  it = ft(:,2);
  ht = ft(:,3);
  ct = ft(:,4);