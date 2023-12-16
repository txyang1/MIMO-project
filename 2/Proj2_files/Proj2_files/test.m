clear all
close all
clc


M = 4 ;
N = 4 ;
H = 1 / sqrt ( 2 )*( randn (M,N) + 1i * randn (M,N) ) ;
P = 10 ^ 3 ;

 % S e t YALMIP o p t i o n s
 options = sdpsettings( 'solver', 'sdpt3', 'verbose', 0);
 % I n i t i a l i z e o p t i m i z a t i o n v a r i a b l e s
 Q = sdpvar(N, N, 'hermitian' , 'complex') ;
 % D e f i n e c o n s t r a i n t s e t
 Constraints = [Q>=0, trace(Q)<=P] ;
 % D e f i n e o b j e c t i v e f u n c t i o n ( m i n i m i z a t i o n o nly ! )
 Objective = - logdet(eye (M) + H*Q*H') ;
 % Solve o p t i m i z a t i o n problem
 sol = optimize(Constraints ,Objective , options) ;
 % R e t r i e v e s o l u t i o n
 C = real(log2(det(eye(M)+ H*value(Q)*H')));