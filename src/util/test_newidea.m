1; % not a function file

function V = projectLow(V)
   V = V(:,1:3);
   V = abs(V);
   for i=1:size(V,2)
	   V(:,i) = V(:,i) - min(V(:,i));
	   V(:,i) = V(:,i) / range(V(:,i));
   end
   V = round(V*100);
end

function [V,v] = getEig(M)
  [V,v] = eig(M);
  [v,idx] = sort(diag(v));
  V = V(:,idx);
end

function L = getLaplacian(M)
  deg = diag(diag(M));
  L = deg - M;
end

function idx = findNN(P)
  idx = zeros(size(P,1),1);
  for i=1:size(P,1)
	  v = zeros(size(P,1),1);
	  for j=1:size(P,1)
		  n = norm(P(i,:)-P(j,:));
		  v(j) = n;
	  end
	  v(i) = inf;
	  [x, ix] = min(v);
	  idx(i) = ix;
  end
end
% --------------------------------------

A = zeros(9,9);

idx  = [ 1 2; 1 4; 2 3; 2 4; 2 6; 3 6; 4 5; 4 8; 5 6; 5 7; 5 8; 6 7; 7 9; 8 9];
%idx = [ 5 6; 5 7; 6 7];

idx2 = [ 5 6; 5 7; 6 7];

for i=1:size(idx,1)
	v = idx(i,:);
	A(v(1), v(2)) = 1;
	A(v(2), v(1)) = 1;
end

X = [];
for i = 1:20
	L = getLaplacian(A);
	[V,v] = getEig(L);
	P = projectLow(V);
	nn = findNN(P);
	X = [X nn];
	for i=1:size(idx2,1)
		v = idx2(i,:);
		A(v(1),v(2)) *= 1.05;
		A(v(2),v(1))  = A(v(1),v(2));
	end
end
