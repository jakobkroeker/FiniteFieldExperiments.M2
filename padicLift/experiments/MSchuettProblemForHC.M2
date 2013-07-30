

restart;

-------------------------- problem construction:

tr := new MutableList ;


#tr#1 = {45,-30,-5,-5,-5 };
#tr#2 = {175,25,-250,-25,75 };
#tr#3 = {750,250,-875,-750,625 };
#tr#4 = {6875,-3125,625,-625,-3750 };
#tr#5 = {-3125,-3125,12500,28125,-34375 };
#tr#6 = {156250,-156250,0,15625,-15625 };

tr#1 = { 45-50, (-30)+50, -5, -5, -5 };
tr#2 = { 175-250, 25, (-250)+250, -25, 75 };
tr#3 = { 750-10*125, 250, -875, -750+10*125, 625 };
tr#4 = { 6875-10*625, -3125, 625, -625, -3750+10*625 };
tr#5 = { -3125, -3125, 12500, 28125, -34375 };
tr#6 = { 156250-10*5^6, -156250+10*5^6, 0, 15625, -15625 };



p := 5;

u := null;
u = symbol u;

varA := flatten toList apply( 7..11, k->
                                        toList apply( 0..p-2, jj->a_{k,jj} ) 
                            );

rng := QQ[varA,u,T]/( u^4 + u^3 + u^2 + u + 1 );

u = null;
u = symbol u;
u = (gens( rng))#20
T = (gens (rng))#21
tl := new MutableList ;
qq := null;


-- all traces @p,p^1,...p^6.


for k in 1..6 do
{
 tl#k=new MutableList;
 for i in 1..p-1 do
 {
    qq = tr#k#0;
    for ii in 1..p-1 do
    {
       qq = qq + tr#k#ii*(u^(i*ii));
    };
    tl#k#i = qq/p;
 };
};

c := new MutableList;
c#0 = new MutableList;
c#1 = new MutableList;

-- Newton-identities

for i from 1 to p-1 do
{
 c#0#i = 1;
 c#1#i = tl#1#i;
};

for k in 2..6 do
{
  c#k = new MutableList;
  for i in 1..p-1 do
  {
    jjIdx := 1..k;
    sumlist = toList apply( jjIdx, jj-> (-1)^(jj-1)*(c#(k-jj)#i)*(tl#jj#i) );
    c#k#i = sum (sumlist)/k;
  };
};


-- remaining coefficients:

for k in 7..11 do
{
  c#k = new MutableList;
  c#k#1 = sum toList apply(p-1, jj-> a_{k,jj}*(u^jj) );
  for i in 2..p-1 do
  {
    c#k#i = sum toList apply( p-1, jj-> (
                             exponent:=sub(mod(i*jj,p),ZZ);
                             if exponent<0 then
                                exponent=exponent+p;
                             a_{k,jj}*(u^exponent)
                         ) );
  };
};

--polynomials:

cp := new MutableList;
for i in 1..4 do
{
 cp#i = sum toList apply (12, jj->(-1)^(jj+1)*T^jj*c#(11-jj)#i );
}

time mm := product( apply (4, i->cp#(i+1) ) );
-- 75 sec.



(M1,C1) := coefficients(mm,Variables=>{T},Monomials=>{T^44})

(M2,C2) := coefficients(mm,Variables=>{T},Monomials=>{T^0});
--print( C1_(0,0)-p^(46-2*1)*C2_(0,0));

idealgens := new List;

-- Relationen 1..22 ( i in {1..22} ) :
-- coeff(mm, T, i-1) = p^(46-2*i)*coeff(  mm, Z, 45-i )

for i in 1..22 do
{
 print i;
 (M1,C1) := coefficients(mm,Variables=>{T},Monomials=>{T^(i-1)});
 (M2,C2) := coefficients(mm,Variables=>{T},Monomials=>{T^(45-i)});
--  print  C1_(0,0)-p^(46-2*i)*C2_(0,0);
 idealgens = append(idealgens, C1_(0,0)-p^(46-2*i)*C2_(0,0) );
}

idgens := ideal idealgens;

rngZ := ZZ[varA,u,T]/ (u^4 + u^3 + u^2 + u + 1);

rngFp := ZZ/2[varA,u,T]/ (u^4 + u^3 + u^2 + u + 1);

idz = clearCoeffDenominators( idgens );

idFp :=  substitute( idz, rngFp );

jac=jacobian idz;

-------------------solutions over F2-----------------------------------------


sol1 := matrix {{  0,0,1,0, 1,1,1,1, 0,1,1,0, 0,0,1,0, 0,0,0,1,0,0 }}
sol1 = sub(sol1, rngFp );
sub( idgens, sol1 ) --- =0
jac:=jacobian idz;
jacsub:=sub(jac,sol1)
rank jacsub --17

sol2 := matrix {{  1,0,0,1, 1,0,1,0, 0,1,0,1, 0,0,1,1, 0,1,1,1 ,0,0}}
sol2 = sub(sol2, rngFp );
jacsub:=sub(jac,sol2)
rank jacsub --20
sub( idgens, sol2)


sol3 := matrix {{  1,0,1,1, 1,1,1,0, 1,1,0,0, 1,0,0,1, 1,0,0,0 ,0,0}}
sol3 = sub(sol3, idgens );
sub( idgens, sol3)
jacsub:=sub(jac,sol3)
rank jacsub --18

sol4 := matrix {{  0,1,1,1, 1,0,0,0, 1,1,0,0, 1,1,0,0, 1,0,0,0, 0,0 }}
sol4 = sub(sol4, rngFp )
sub( idgens, sol4)

jac:=jacobian idgens;


jacsub:=sub(jac,sol4)
rank jacsub --20

sol5 := matrix {{  1,0,1,1, 0,1,1,0, 0,0,1,0, 1,1,0,1, 1,0,1,0, 0,0 }}
sol5 = sub(sol5, rngFp );
sub( idgens, sol5)
jacsub:=sub(jac,sol5)
rank jacsub --17
sol6 := matrix {{  0,1,0,1, 0,0,0,0, 0,1,1,0, 0,0,0,1, 1,0,1,1, 0,0 }}
sol6 = sub(sol6, rngFp );
sub( idgens, sol6)
jacsub:=sub(jac,sol6)
rank jacsub --15


sol7 := matrix {{  1,0,1,1, 0,1,1,0, 0,1,0,1, 1,1,1,0, 1,0,1,1 , 0,0 }}
sol7 = sub(sol7, rngFp );
sub( idgens, sol7)
jacsub:=sub(jac,sol7)
rank jacsub --17


sol8 := matrix {{  0,0,0,0, 1,1,0,0, 1,0,1,0, 0,0,0,1, 1,1,1,1 ,0,0}}
sol8 = sub(sol8, rngFp );
sub( idgens, sol8)
jacsub:=sub(jac,sol8)
rank jacsub --12


----------------------------- prepare lifting :

loadPackage ("padicLift",Reload=>true)
installPackage("padicLift")

--varAZ:=(gens rngZ)_{0..19} --skip variables 'T' and 'u'

varAZ:=(gens ring idz)_{0..19} --skip variables 'T' and 'u'
rngred :=  ZZ[varAZ]
idzred := sub (idz,rngred); -- problem ideal in ZZ 

-- Lösungspunkte in  rngred 

sol2red := matrix {{  1,0,0,1, 1,0,1,0, 0,1,0,1, 0,0,1,1, 0,1,1,1 }}
sol2red = sub(sol2red, rngFp )
sol4red := matrix {{  0,1,1,1, 1,0,0,0, 1,1,0,0, 1,1,0,0, 1,0,0,0 }}
sol4red = sub(sol4red, rngFp )

 
solutions={matrix {{0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0,
0, 0}}, matrix {{1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
0}}, matrix {{0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0,
0}}, matrix {{1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 0,
0}}, matrix {{1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1,
0}}, matrix {{0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0,
1}}, matrix {{1, 1, 1, 1, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1,
1}}, matrix {{0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
1}}}
solutions=apply(solutions, sol->sub(sol,ZZ/2))

sub( idzred, sol2red ) -- =0
sub( idzred, sol4red ) -- =0

varAZred=(gens rngred)_{0..19}

-- jac = jacobian idzred;
-- jacsubred = sub(jac, sol4red)

ranks= apply(solutions, sol-> rank sub(jac, sol))


--------------------------- lift: 
lllres := computeMinPolys( idzred, solutions#6, varAZred);

lllres := computeMinPolys( idzred, sol2red, varAZred);

lllres := computeMinPolys( idzred, sol4red, varAZred);





