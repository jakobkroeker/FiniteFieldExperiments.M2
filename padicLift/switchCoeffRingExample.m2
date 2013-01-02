--substituteExample

-- test: what happens, if we use change the coefficient ring 
--       to ZZ mod p^k ( some problems are caused by ZZ[]/(prime^k) -rings !)


 loadPackage("padicLift",Reload=>true);
   installPackage("padicLift")
  
  x:=null; x=symbol x;

   prime := 7;

    RQ := QQ[x];
    x := (gens(RQ))#0;
    FQ := (1/3*x+1/2) ; 
  
    IFQ := ideal FQ;
    IFZ := disposeRationalCoeffs(IFQ);  

    Fp := ZZ/prime[x];
    Fpe := ZZ[]/prime[x];

    RZ := ZZ[x];
    
    -- substitute(IFQ,RZ) does not what we probably would expect ( ideal(1/3*x+1) over QQ <=> ideal(x+3) over ZZ )

     substitute(IFQ,RZ);


    --cannot promote 1/3 to Fp, but substitute does not care.

    try ( promote(1/3,Fp)) then (error("this should fail") )  else ( );

    substitute(IFQ,Fp);

     --disposeRationalCoeffs does the job correctly, but only if we have not a nested ring or a quotient ring. 
    IFZ := disposeRationalCoeffs(IFQ);

   IFP := sub(IFQ, Fp); --correct
  IFPe := sub(IFQ, Fpe); --wrong

   IFP := sub(IFZ, Fp); --correct
  IFPe := sub(IFZ, Fpe);--correct

   solutionQQ := matrix{{-3/2}} ;
   solution := sub( solutionQQ ,ZZ/7) ;
   sub( IFQ, solution);
   sub( IFZ, solution  );

   nextApprox := nextLift(IFQ,solution);
   sub(IFQ,nextApprox);
  sub(IFZ,nextApprox);

 solutionQQ := matrix{{-3/2}} ;


  apply(prime-1, i->(
   nextApprox:=sub(  matrix{{i+1}}  ,ZZ[]/prime ) ;
    ( sub(IFQ,nextApprox), sub(IFZ,nextApprox), sub(IFP,nextApprox), sub(IFPe,nextApprox) )
))

  apply(prime-1, i->(
   nextApprox:=sub(  matrix{{i+1}}  ,ZZ/prime ) ;
    ( sub(IFQ,nextApprox) , sub(IFZ,nextApprox) , sub(IFP,nextApprox), sub(IFPe,nextApprox) )
))

  apply((prime^2-1, i->(
   nextApprox:=sub(  matrix{{i+1}}  ,ZZ[]/prime^2) ;
    sub(IFQ,nextApprox)
))

 

    x := (gens ring IFQ)_0;

   -- wrong result!

   result := computeMinPolys(IFQ, solution,{x} );

    x := (gens ring IFZ)_0;

   correctResult := computeMinPolys(IFZ, solution,{x} );


