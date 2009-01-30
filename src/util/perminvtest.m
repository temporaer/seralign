% vim:ft=matlab
function [A,B,va,vb,M] = perminvtest(n)
  kernelp = 2;
  M = get_rand_adj2(n);

  La = get_laplacian(M);
  Ka = getKernel(La, kernelp);
  [A,va] = seig(La);

  P = get_perm(M);
  M = P*M*P';
  Lb = get_laplacian(M);
  Kb = getKernel(Lb, kernelp);
  [B,vb] = seig(Lb);

  Ha = heatkernel(A,va, 0.0001);
  Hb = heatkernel(B,vb, 0.0001);
  %A=Ha;B=Hb;
  A=round(A*100);
  B=round(B*100);
end

function [V,v] = seig(M)
  [V,v] = eig(M);

  % sort cols according to eigenvalues
  [v,idx] = sort(diag(v));
  V = V(:,idx);

  % sort lines according to fiedler vec
  V = abs(V);
  %[s,idx]=sort(V(:,2));
  %V = V(idx,:);
end

function M = get_rand_adj(n)
  M = rand(n,n) - 0.5;    % make 0/1 random matrix
  M(M<0) = 0;
  M(M>0) = 1;

  for i=1:n               % mirror triangles
    for j=i:n
	  M(i,j) = M(j,i);
	end
  end

  M = M - diag(diag(M));  % make 0s on diagonal
end

function M = get_rand_adj2(n)
  M = [
  0   1   1   0   0   0;
  1   0   1   0   0   0;
  1   1   0   .1   0   0;
  0   0   .1   0   1   1;
  0   0   0   1   0   1
  0   0   0   1   1   0
  ];
end

function L = get_laplacian(M)
  n = size(M,1);
  D = diag(sum(M)); % degree matrix
  L = D - M;        % laplacian
  L = eye(n,n) - D^-.5 * M * D^-.5;   % normalized laplacian
  %L = M;  % adjmat
end

function K = getKernel(L, p)
  % p-step random walk
  n = size(L,1);
  K = (2*eye(n,n) - L)^p;
end

function P = get_perm(M)
  n = size(M,1);
  r = rand(n);
  [s, idx] = sort(r);
  P = zeros(n,n);
  for i=1:n
    P(i,idx(i)) = 1;
  end
end
