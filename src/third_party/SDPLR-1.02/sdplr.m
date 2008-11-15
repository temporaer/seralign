%
%  [x,y,info,r] = sdplr(A,b,c,K,pars,lrA)
%
%  SDPLR 1.02      http://dollar.biz.uiowa.edu/~sburer/software/SDPLR
%                  Email and bug reports to: samuel-burer@uiowa.edu
%
%  > x = sdplr(A,b,c,K)
%
%    Solves the semidefinite program:   min  c'*x  st  A*x = b, x in K
%
%    K   describes  the   cone  of   nonnegative  linear   and  positive
%    semidefinite variables. K follows the SeDuMi format (see the SeDuMi
%    website http://sedumi.mcmaster.ca/).
%
%    K  can  have two  fields,  K.l  and  K.s, representing  Linear  and
%    Semidefinite. K.l is the number of nonnegative components of x, and
%    K.s lists the dimensions of the semidefinite constraints on x.
%
%    Example: K.l = 10,  K.s = [4 3], and K.l listed  before K.s. Then x
%    has dimension 35 = 10 + 4*4 + 3*3 and satisfies the constraints
%
%                   x(1:10)                nonnegative
%                   reshape(x(11:26),4,4)  positive semidefinite
%                   reshape(x(27:35),3,3)  positive semidefinite
%
%    Example: K.s = [3  4], K.l = 10, and K.s listed  before K.l. Then x
%    has dimension 35 = 10 + 4*4 + 3*3 and satisfies the constraints
%
%                   reshape(x(1:9),3,3)    positive semidefinite
%                   reshape(x(10:25),4,4)  positive semidefinite
%                   x(26:35)               nonnegative
%
%  > [x,y] = sdplr(A,b,c,K)
%
%    Returns the dual solution of:   max  b'*y  st  c - A'*y in K
%
%  > [x,y,info] = sdplr(A,b,c,K)
%
%    Returns information from the execution of SDPLR, where:
%
%        info.majoriter  =  # of major iterations (i.e., penalty increases)
%        info.minoriter  =  # of minor iterations (i.e., steps)
%        info.cgs        =  # of conjugate gradient iterations
%        info.time       =  total time (in seconds)
%        info.penalty    =  final value of penalty parameter
%
%  > [x,y,info,r] = sdplr(A,b,c,K)
%
%    Returns the internal format r of x used by SDPLR.
%
%    The variable r  is a cell array, such that  r{i} stores the portion
%    of x correpsonding to the i-th cone constraint on x, as given by K:
%
%      portion x  =  r{i}.*r{i}                    if linear
%      portion x  =  reshape(r{i}*r{i}',sz*sz,1)   if semidef of size sz
%
%    Example: K.l = 10, K.s = [4 3], and K.l listed before K.s. Then
%
%                   x(1:10)   =  r{1}.*r{1}
%                   x(11:26)  =  reshape(r{2}*r{2}',9,1)
%                   x(27:35)  =  reshape(r{3}*r{3}',16,1)
%
%  > [x,y,info,r] = sdplr(A,b,c,K,pars)
%
%    Passes paramaters into SDPLR, where:
%
%        pars.feastol     =  Desired level of infeasibility in A*x = b
%                            (default = 1.0e-5)
%
%        pars.centol      =  Desired level of accuracy in intermediate
%                            calculations (default = 1.0e-1)
%
%        pars.dir         =  Method for calculating step direction
%                            (1 = limited memory BFGS, 2 = truncated Newton)
%                            (default = 1)
%
%        pars.penfac      =  Factor by which the penalty parameter is
%                            increased each major iteration (default = 2.0)
%
%        pars.reduce      =  Whether or not to reduce the problem dimension-
%                            ality as the algorithm progresses (1 = yes,
%                            0 = no) (default = 1)
%
%        pars.limit       =  Limit on the total time (in seconds)
%                            (default = 3600)
%
%        pars.printlevel  =  Whether to print output (0 = no, 1 = yes)
%                            (default = 1)
%
%   Other,  less  critical  parameters  are  available.  Please  contact
%   samuel-burer@uiowa.edu for more information.
%
%   To save on memory, there is also:
%
%        pars.soln_factored  =  Return x in its factored form (i.e., x <- r)
%                               (1 = yes, 0 = no) (default = 0)
%
%  > [x,y,info,r] = sdplr(A,b,c,K,pars,lrA)
%
%    Passes low-rank data matrices into SDPLR,  where lrA is an array of
%    structures.  E.g.,  lrA =  [lrA(1),  lrA(2),  lrA(3)] passes  three
%    low-rank matices.
%
%    Each low-rank  matrix must apply to  a (semidefinite) cone i  and a
%    constraint j,  where the cones are  given by K and  the constraints
%    correspond to  the rows of A.  One can also indicate  the objective
%    function by j=0. The matrix is  given by V*diag(D)*V', where V is a
%    matrix and D  is a vector. The number  of rows of V is  the size of
%    the semidefinite cone, and the number of columns of V (also size of
%    D) is the rank of the matrix.
%
%    Each structure in lrA has four fields:
%
%        cons   =  which constraint j
%        start  =  index of the first position of cone i in x
%        D      =  vector D
%        V      =  matrix V
%
%    Example: K.l = 10, K.s = [4 3],  and K.l listed before K.s. If V is
%    size (3,2) and D is size 2, then the rank-2 matrix V*diag(D)*V' can
%    be applied to  x in the third cone (semidef)  and second constraint
%    as
%
%         lrA(1).cons   =  2 
%         lrA(1).start  =  27
%         lrA(1).D      =  D
%         lrA(1).V      =  V
%
%    Note: Even if the j-th constraint is specified completely by
%          lrA, it is still necessary for A to have a j-th row (sparse
%          and empty) to serve as a placeholder.
%
%  ------------------------------------------------------------------------
%
%  SDPLR 1.02      http://dollar.biz.uiowa.edu/~sburer/software/SDPLR
%                  Email and bug reports to: samuel-burer@uiowa.edu
function [varargout] = sdplr(varargin);

if nargout == 0
    mexsdplr(varargin{:});
else
    [varargout{1:nargout}] = mexsdplr(varargin{:});
end
