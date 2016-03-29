% vardec.m: Uses the parameter estimates and standard errors 
%             found by estse.m to compute variance decompositions
%             and standard errors for the real business cycle
%             model with indivisible labor.
%
%           Returns a 56x2 matrix vdec.  The first column is:
%
%             [ yvar cvar ivar hvar yvart cvart ivart hvart ]'
%
%             where yvar, cvar, ivar, and hvar are vectors of k-step
%             ahead forecast error variances in output, consumption,
%             investment, and hours worked, and k = 1,4,8,12,20,40,
%             and infinity.  The variances are all expressed as 
%             percentages.  The vectors yvart, cvart, ivart, and
%             hvart are vectors of percentages due to the technology
%             shock.  The second column gives the standard errors
%             for the corresponding elements of the first column.
%
%           Formulas used to compute variance decompositions and
%             standard errors are contained in vardecfn.m.
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

% find covariance matrix of estimated parameters

  bighx = inv(hessmat);

% compute variance decompositions

  vdec1 = vardecfn(tstar);

% compute standard errors

  epsmat = 1e-6*eye(21);

  gradmat = zeros(56,21);

  for i = 1:21

    gradmat(:,i) = (vardecfn(tstar+epsmat(:,i))-vdec1)/1e-6;

  end

  vdec2 = gradmat*bighx*gradmat';

  vdec2 = sqrt(diag(vdec2));

% collect results

  vdec = [ vdec1 vdec2 ];