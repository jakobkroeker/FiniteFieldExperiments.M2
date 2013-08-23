
tr := new MutableList ;

tl := new MutableList ;

tr#1 = {45,-30,-5,-5,-5};
tr#2 = {175,25,-250,-25,75};
tr#3 = {750,250,-875,-750,625};
tr#4 = {6875,-3125,625,-625,-3750};
tr#5 = {-3125,-3125,12500,28125,-34375};
tr#6 = {156250,-156250,0,15625,-15625};

tr#1 = { 45-50, (-30)+50, -5, -5, -5 };
tr#2 = { 175-250, 25, (-250)+250, -25, 75 };
tr#3 = { 750-10*125, 250, -875, -750+10*125, 625 };
tr#4 = { 6875-10*625, -3125, 625, -625, -3750+10*625 };
tr#5 = { -3125, -3125, 12500, 28125, -34375 };
tr#6 = { 156250-10*5^6, -156250+10*5^6, 0, 15625, -15625 };


p := 5;

u := null;
u := symbol u;

varA := flatten toList apply(7..11, k->toList apply(0..p-2, jj->a_{k,jj}));

rng2 := QQ[u,T,varA];



u := (gens( rng2))#0;

idu := ideal (u^4 + u^3 + u^2 + u + 1);


for k in 1..6 do { tl#k=new MutableList; for i in 1..p-1 do {
  qq=tr#k#0;   for ii in 1..p-1 do  {     qq=qq + tr#k#ii*(idu^(i*ii));    };
  tl#k#i=qq/p; }; };

c := new MutableList;
c#0 = new MutableList;
c#1 = new MutableList;

for i from 1 to p-1 do {  c#0#i=1;  c#1#i=tl#1#i; };

for k in 2..6 do  { c#k = new MutableList; for i in 1..p-1 do
{ jjIdx := 1..k; sumlist = toList apply( jjIdx, jj-> (-1)^(jj-1)*(c#(k-jj)#i)*(tl#jj#i) );
c#k#i = sum (sumlist)/k; }; };





for k in 7..11 do  {  c#k = new MutableList;
  c#k#1 = sum toList apply(p-1, jj-> a_{k,jj}*(idu^jj) );
  for i in 2..p-1 do   {   c#k#i = sum toList apply( p-1, jj-> (
      exponent := sub(mod(i*jj,p),ZZ);       if exponent<0 then      exponent = exponent+p;
      a_{k,jj}*(idu^exponent)     )     );  }; };

cp := new MutableList;
for i in 1..4 do { cp#i = sum toList apply (12, jj->(-1)^(jj+1)*T^jj*c#(11-jj)#i ); }


mm := product( apply (4, i->cp#(i+1) ) );


(M1,C1) := coefficients(mm,Variables=>{T},Monomials=>{T^44});

(M2,C2) := coefficients(mm,Variables=>{T},Monomials=>{T^0});

idealgens:= new List;

for i in 1..22 do {  print i;
 (M1,C1) := coefficients(mm,Variables=>{T},Monomials=>{T^(i-1)});
 (M2,C2) := coefficients(mm,Variables=>{T},Monomials=>{T^(45-i)});
 idealgens = append(idealgens, C1_(0,0)-p^(46-2*i)*C2_(0,0) ); }


idealgens = append( idealgens, u^4 + u^3 + u^2 + u + 1);

id:= ideal idealgens;

rngZ := ZZ[u,T,varA];

rngFp := ZZ/2[u,T,varA];


idz := substitute(id,rngZ);

idFp :=  substitute(id,rngFp);

jac := jacobian idFp;
--sol:= matrix {{ 1,0, 0,0,1,0, 1,1,1,1, 0,1,1,0, 0,0,1,0, 0,0,0,1 }}

--#sol=sub(sol,ZZ/2);

--substituted:=sub(id,sol);
--sub(substituted,ZZ/2);

--sol:= matrix {{ 0,0, 0,1,0,0, 1,1,1,1, 0,1,1,0, 0,1,0,0, 1,0,0,0 }}

--sol=sub(sol,ZZ/2);

--sub(id,sol);

--sol:= matrix {{ 0,0, 0,0,1,0, 1,1,1,1, 0,1,1,0, 0,0,1,0, 0,0,0,1 }}

--sol=sub(sol,ZZ/2);

--sub(id,sol);

--zerosol:= matrix {{ 0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 }}

--zerosol=sub(sol,ZZ/2);

--sub(sol,id);

rngZ := ZZ[u,T,varA];

rngFp := ZZ/2[u,T,varA];


idz := substitute(id,rngZ);

idFp :=  substitute(id,rngFp);

jac := jacobian idFp;

rank jac

