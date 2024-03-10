using TT

d = 60;
n = 2^d;
h = 1.0/(n-1)

x = qtt_x( 2 .* ones(Int,d) );
tt = funcrs(x,sqrt,1e-12,x,8);
p = tt_ones( 2 .* ones(Int,d) );

error = dot( p,tt )*h - 2/3 # Difference between numerical and theoretica