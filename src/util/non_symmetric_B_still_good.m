#!/usr/bin/octave


A = rand(10,10);
%A = A>0.5;
O = 2*eye(size(A,1)); O(1) = 1; O(end) = 1;
Op = O^0.5;
Om = O^-0.5;


for i=1:size(A,1)
	for j=i:size(A,2)
		A(j,i) = A(i,j);
	end
end

org = Op * A * Om;
sym = Op * A * Om;

for i=1:size(A,1)
	for j=i:size(A,2)
		sym(j,i) = sym(i,j);
	end
end

orgval = [];
symval = [];

fact = 1;
while(length(orgval)< 50)
	R = fact*rand(size(A,1), size(A,2))-(fact/2);
	isposdef = true;
	for i=1:length(R)
		if ( det( R(1:i, 1:i) ) <= 0 )
			isposdef = false;
			break;
		end
	end
	if (isposdef)
		printf('.');fflush(stdout);
		orgval = [ orgval trace(org * R) ];
		symval = [ symval trace(sym * R) ];
	end
end

hold off;
plot(orgval,symval,'+');
hold on;
plot([min([orgval symval]) max([orgval symval])], [min([symval orgval]), max([symval orgval])])
xlabel('trace (B * Y)')
ylabel('trace (B'' * Y)')
