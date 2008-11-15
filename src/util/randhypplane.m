#!/usr/bin/octave
# vim:ft=matlab

% sdp-soln
values
n = size(A,1);

% functions
Bs = [];
function Bs = store_b (Bs,B)
	Bs = [Bs ; B ];
end
function B = get_b (Bs,i,n)
	B = Bs((i-1)*n+1:(i*n)-1, :);
end

% omega
O = 2*eye(n); O(1) = 1; O(end) = 1;
Op = O^0.5;
Om = O^-0.5;


org = Op * A * Om;
sym = Op * A * Om;

for i=1:size(A,1)
	for j=i:size(A,2)
		sym(j,i) = sym(i,j);
	end
end


if 1 == 1

	X  = Y;
	Xp = Y;  % X = X + sum ... -- pertubed.

	for i=1:100
		R  = chol(Xp);
		rho = rank(Xp);

		% Z = in S^rho    s.t. <Z, R' Aj R> = 0    forall j=1..m. Here: Aj=eye(n)
		Z = R' * eye(n) * R * rand(n,n);

		% phi(Z) = rank(Z)==rho ? 1 : -1
		L = eig(Z);
		if sum(diag(L)<-0.000000001) > 0
			phiZ = -1;
		else
			phiZ = 1;
		end
		
		ti_inv = max( phiZ * L )
		
		if ti_inv == 0
			break
		end

	end


	if 0 ==  1

		y = ys(:,i);
		x = Op \ y;

		done = zeros(n,1);
		path = [];

		x = x + min(x) + 1;
		[m,i] = max(x);
		while(size(path) < n)
			done(i) = 1;
			path = [path i];
			x(i) = -1;
			j    = i;
			[m,i] = max(x.*A(:,i));
			if( (done(i) || A(i,j) <1 ) && size(path) < n)
				[m,i] = max(x);
				printf('Jumping from %d to %d...\n', j,i);
			end
		end
	end
end
