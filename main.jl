using TT
A = randn(5,5,5,5,5);
Att = VectorTT(A,eps());
A = randn(5,5,5,5,5,5);
Attm = MatrixTT(A);
Qttm = randn(Att);