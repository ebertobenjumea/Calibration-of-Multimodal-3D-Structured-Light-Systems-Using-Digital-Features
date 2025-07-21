x = pout(1, :); y = pout(2,:); X = pin(1,:); Y = pin(2,:);
rows0 = zeros(3, n);
rowsXY = -[X; Y; ones(1,n)];
hx = [rowsXY; rows0; x.*X; x.*Y; x];
hy = [rows0; rowsXY; y.*X; y.*Y; y];
h = [hx hy];


for i=1:n
    X=pin(1,i); Y=pin(2,i); %Z=pin(3,i);
    x=pout(1,i); y=pout(2,i);
    A(2*(i-1)+1,:)=[-X -Y -1 0 0 0 x*X x*Y x];
    A(2*(i-1)+2,:)=[0 0 0 -X -Y -1 y*X y*Y y ];
end
[U,D,V] = svd(A);
x=V(:,end) 
% Reescribimos x en forma de la matriz M
M=[ x(1) x(2) x(3);
x(4) x(5) x(6);
x(7) x(8) x(9)];