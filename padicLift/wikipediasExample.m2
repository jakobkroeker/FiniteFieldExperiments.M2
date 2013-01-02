
-- Hensel Lifting example from http://en.wikipedia.org/wiki/Hensel%27s_lemma
x := null;
x = symbol x;
rng := QQ[x];

pol := x^2+1;
IP := ideal pol;
prime := 2;
Fp := ZZ/prime[x];

Fp4 := ZZ[]/prime^2[x];

point := matrix {{1_Fp}} ;
sub(IP,point)
assert ( sub(IP,point)==0 );

JZ := jacobian IP

JFp := sub(JZ,point) -- is zero

JZS := sub(JZ,sub(point,ZZ))

sub(IP, sub(point,Fp4 )) --is not zero => there is no lifting!

-----------------------------------------------------------------------------

-- Hensel Lifting Example from http://en.wikipedia.org/wiki/Hensel%27s_lemma
x := null;
x = symbol x;
rng := QQ[x];

pol := x^2-17;
IP := ideal pol;

prime := 2;
Fp := ZZ/prime[x];
Fp4 := ZZ[]/prime^2[x];

point := matrix {{1_Fp}} ;
sub(IP,point)
assert ( sub(IP,point)==0 );
assert ( sub(IP, sub(point,Fp4 ) )==0 );

JZ := jacobian IP

JFp := sub(JZ,point)



sub( IPZ, sub( point,Fp4 ) )==0 

IPZ := disposeRationalCoeffs( IP); 
assert ( sub(IPZ,point)==0 );
assert ( sub(IPZ, sub(point,Fp4 ) )==0 );
 
liftOptions= new LiftOptions
liftOptions#"minColumnNormDistanceFactor"=10
liftOptions#"verbose"=true
computeMinPolys (IPZ,point,{x},"liftAndLLLOptions"=>liftOptions)

point = nextLift(IPZ,point)


-- jacobian is zero and F(point) mod p^(k+1) is zero => every point+ t*prime^k (k=1) is a solution:
-- this is not implemented, since I do currently not know, how to extrapolate a multidimensional version.


nextSolution :=  random(ZZ^1,ZZ^1)*matrix{{ prime }} + sub( point,Fp4 );
sub( IP, nextSolution )
assert( sub( IP, nextSolution )==0 );

1*matrix{{2}}+sub(point,Fp4 ) 

-- approach from Anton Leykin:

r = rank JFp
-- r=0  =>r+1=1: 

h = random (ZZ^1,ZZ^1); -- a (   r+1 )-vector
B = random (ZZ^1,ZZ^1);  -- a ( n × r+1 )-matrix

C = JZ*B;
CP := sub( C ,point);

while (rank CP!=0) do (
h = random (ZZ^1,ZZ^1); 
B = random (ZZ^1,ZZ^1);
C = JZ*B;
CP := sub( C ,point);
);

CP := sub( C ,sub(point,QQ) );

newrng := QQ[x,lambda];
x = (gens newrng)#0;
lambda = (gens newrng)#1;
pol := x^2-17;
neweq1= lambda*sub(CP_(0,0),newrng)
neweq2= h_(0,0)*lambda-1


NIP = ideal {pol, neweq1, neweq2 }

NJ = jacobian NIP

npoint := matrix {{1_Fp,1_Fp}}

sub(NIP,npoint)

  NIPZ := disposeRationalCoeffs(NIP); 
sub(NIPZ,npoint)

nextLift( NIPZ, npoint);





