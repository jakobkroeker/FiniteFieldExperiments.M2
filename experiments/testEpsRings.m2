-- test performance of eps-rings

restart
K = ZZ/11
E = K[e]
E2 = K[f]/f^2
R = K[x,y,z,w]

pointE = random(K^1,K^4)+e*sub(random(K^1,K^4),E)
pointF = sub(pointE,e=>f)

m=10
M = random(R^{m:0},R^{m:-1})
time (detM = det M);
-- used 0.823192 seconds (m=10)
time det sub(M,pointE)
-- used 0.001149 seconds (m=10)
time det sub(M,pointF)
-- used 7.10316 seconds (!!!!)
time det(sub(M,pointF),Strategy=>Bareiss)
-- used 0.003823 seconds
time det(sub(M,pointF),Strategy=>Cofactor)
 -- used 5.73505 seconds
 
time sub(detM,pointE) 
-- used 0.003164 seconds (m=10)
time sub(detM,pointF)
-- used 0.002314 seconds (m=10)


restart
K = ZZ/11
E = K[e]
d = 100
E2 = K[f]/f^d
R = K[x,y,z,w]

pointE = sum apply(d,i->e^i*sub(random(K^1,K^4),E));
pointF = sub(pointE,e=>f);
point = sub(pointE,K)

m=10
M = random(R^{m:0},R^{m:-1})
time (detM = det M);
-- used 0.823192 seconds 
time det sub(M,pointE);
-- used 22.7502 seconds (d=100)
time det(sub(M,pointF),Strategy=>Bareiss);
-- used 2.32395 seconds (d=100)

time sub(detM,pointE);
-- used 26.6599 seconds (d=100)
time sub(detM,pointF);
-- used 4.13755 seconds (d=100)
time sub(detM,point)
-- used 0.00023 seconds

m=10
M = random(R^{m:0},R^{m:-2})
time (detM = det M);
-- used 27.9315 seconds
time det sub(M,pointE);
-- used 97.9807 seconds (d=100)
time det(sub(M,pointE),Strategy=>Bareiss);
-- used 98.277 seconds
time det(sub(M,pointF),Strategy=>Bareiss);
-- used 4.4757 seconds (d=100)

betti basis(10,R)
apropos "LR*"
help LUdecomposition
points = apply(20,i->random(K^1,K^4))
points0 = select(points,point->0==det sub(M,point))
point0 = points0#0
LUdecomposition (sub(M,point0))
LU = LUdecomposition (sub(M,point))
LU#1^-1
LU#2^-1
apropos "code"
code LUdecomposition
help LUdecomposition

LU = LUdecomposition(random(K^10,K^10)*diagonalMatrix{1,1,1,1,1,1,1,1,0,0}*random(K^10,K^10))
LU#1*LU#2


--- evaluation
restart
K = ZZ/11
E = K[e]
d = 2
E2 = K[f]/f^d
R = K[x,y,z,w]

pointE = sum apply(d,i->e^i*sub(random(K^1,K^4),E));
pointF = sub(pointE,e=>f);
point = sub(pointE,K)


time F = random(60,R);
-- used 1.84935 seconds (30)
-- used 14.0919 seconds (40)
-- used 76.9105 seconds (50)
-- used 239.726 seconds (60)
--time sub(F,pointE);
time sub(F,pointF)
-- used 0.069715 seconds (30)
-- used 0.17217  seconds (40)
-- used 0.407939 seconds (50)
-- used 0.789915 seconds (60)
time sub(F,point)
-- used 0.000215 seconds (30)
-- used 0.000244 seconds (40)
-- used 0.000789 seconds (50)
-- used 0.000306 seconds (60)


time F = random(60,R);
trials = 10000
time (sum apply(trials,i->(
	  point = random(K^1,K^4);
	  (timing sub(F,point))#0
	  )))/trials
-- 0.00025 (40)
-- 0.00045 (50)
-- 0.00095 (60)
d=2
trials = 100
(sum apply(trials,i->(
	  pointF =  sum apply(d,i->f^i*sub(random(K^1,K^4),E2));
	  (timing sub(F,pointF))#0
	  )))/trials

-- 0.23  (40)
-- 0.44  (50)
-- 0.83  (60)
