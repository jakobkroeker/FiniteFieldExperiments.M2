newPackage(
     "BlackBoxIdeals",
     Version => "0.1", 
     Date => "15.02.2013",
     Authors => {{
           Name => "Jakob Kroeker", 
           Email => "kroeker@uni-math.gwdg.de", 
           HomePage => "http://www.crcg.de/wiki/User:Kroeker"},{
           Name => "Hans-Christian Graf v. Bothmer", 
           Email => "bothmer@math.uni-hannover.de", 
           HomePage => "http://www.crcg.de/wiki/Bothmer"}    
      },
     Configuration => {},
     PackageExports => {"M2Logging"},
     Headline => "black boxes for explicit and implicitly given ideals",
     DebuggingMode => true
)


   -- Q: how to get the parameter list of a function or method?
   --
   -- answer:  not possible yet, but it is possible to get the number of parameters 
   -- using 'disassemble' for functions (where the number of variables is fixed and not variable)
   -- and for methods it forks as follows: consider a method foo with several functions installed.
   -- then   apply(methods foo, m-> (m, disassemble lookup m) )



needsPackage "M2Logging";


export {
    --clearCoeffDenominators,
    "clearCoeffDenominators",
    BlackBoxIdeal,
    BlackBoxParameterSpace,
    blackBoxParameterSpace,
    blackBoxIdeal,
    blackBoxIdealFromEvaluation,
    BlackBoxLogger,
    getEpsRing,
    jetAt,
    isCertainlySingularAt,
    isProbablySmoothAt,
    "keysWithoutSymbols",
    "rebuildBlackBox",
    "guessAcceptedParameterNumber",
    "dropDegreeInfo"
}

--undocumented {
--
--}

-- why is/was setting 'jacobianAt' necessary??? (jk)

-- jacobianAt := global jacobianAt;
-- jacobianAt = global jacobianAt;

idealBlackBoxesProtect = ()->
(
--protect eps;
protect jacobianAt;
protect rankJacobianAt;
protect valuesAt;
protect unknownIsValid;
protect numVariables;
protect numGenerators;
protect isZeroAt;
protect pointProperties;
protect registerPointProperty;
protect rpp;
protect numTrials;
protect setSingularityTestOptions;
 
protect setPointProperty;
protect setValuesAt;
protect checkInputPoint;
protect deduceNumGenerators;
protect setIsZeroAt;
protect getUpdatedBlackBox;
protect type;
protect unknowns;
protect pointProperty;
protect updatePointProperty;
protect setJacobianAt;
protect equations;

protect hasPointProperty;
protect knownProperties;
protect knownPointProperties;
protect knownPointPropertiesAsSymbols;
protect knownMethods;
protect knownAttributes;
);

--todo: fix dublicate code,  -  padicLiftProtect and padicLiftExport

idealBlackBoxesExport = ()->
(
    exportMutable( eps );
    exportMutable( jacobianAt);
    exportMutable( rankJacobianAt);
    exportMutable( valuesAt);
    exportMutable( unknownIsValid);                 
    exportMutable( numVariables);  
    exportMutable( numGenerators);
    exportMutable( isZeroAt);      
    exportMutable( pointProperties);  
    exportMutable( registerPointProperty); 
    exportMutable(rpp);
   exportMutable(upp);
    exportMutable( setPointProperty);
    exportMutable( setValuesAt);    
    exportMutable(  checkInputPoint);
    exportMutable( deduceNumGenerators );
    exportMutable( setIsZeroAt );
   exportMutable( getUpdatedBlackBox );
   exportMutable( type );
 exportMutable( unknowns );
 exportMutable( pointProperty );
 exportMutable( updatePointProperty );
 exportMutable( setJacobianAt );
 exportMutable( equations );

 exportMutable( hasPointProperty );
  exportMutable( knownPointProperties );
  exportMutable (knownProperties);
  exportMutable( knownPointPropertiesAsSymbols );

  exportMutable(knownMethods);
 exportMutable(knownAttributes);
  exportMutable(numTrials);
 exportMutable(setSingularityTestOptions);

)

-- the following symbols which are marked as undocumented are in fact documented inside the BlackBoxParameterSpace and BlackBoxIdeal
-- please do only mark documented symbols as undocumented, at least there should be a comment note inside this package.
-- 
undocumented { 
    setIsZeroAt,      -- internal   
    setJacobianAt,    -- internal   
    setValuesAt,      -- internal   
    setPointProperty, -- internal   
    deduceNumGenerators, -- internal   
    dropDegreeInfo,      -- internal   
    getUpdatedBlackBox,
    pointProperties,     -- internal key     
    equations,
    keysWithoutSymbols,
    checkInputPoint,
    knownPointPropertiesAsSymbols,
    type,
    unknownIsValid,
    unknowns,
    numTrials,
    pointProperty -- internal function
} 


-- swith between protect and export - both are not possible!

--idealBlackBoxesProtect() -- protect the symbols for checking package correctness: no symbol variables should be overwritten by accident!
idealBlackBoxesExport(); -- export the symbols to make the package work 


needsPackage "SimpleDoc";
needsPackage "Text";


BlackBoxLogger = Logger("BlackBoxIdeals");

-- todo: how to switch this on and off?
-- if BlackBoxIdeals#Options#DebuggingMode then 
--    BlackBoxLogger.setLogLevel(LogLevel.DEBUG);

if BlackBoxIdeals#Options#DebuggingMode then
  errorDepth=0;



bblog := BlackBoxLogger;

-- why did I need comparison between a string and a symbol ? (jk)

String ? Symbol := (str,symb)->
(
  1 ? 2
);

Symbol ? String := (str,symb)->
(
  2 ? 1
);


-- find out for a function, how many parameters it does accept.
-- if a function accepts variable number of parameters, returns null
-- if it did not find 'numparms:' in the disasseble string, returns null.
guessAcceptedParameterNumber = method();
guessAcceptedParameterNumber( Function ) := ZZ => (foo)->
(
    lst := disassemble foo;

    bblog.debug  ("disassemble result: " | lst );
    lst = separate( " ", lst );

     restargsPos := position( lst, (str)-> str=="restargs:" );
     numparmsPos := position( lst, (str)-> str=="numparms:" );



     if restargsPos=!=null then 
     (
         newlst := drop (lst, restargsPos);
         restargsPos = position( newlst, (str)-> str=="restargs:" ); 
         
         if restargsPos=!=null then 
         (
           bblog.info ("do not know how to handle methods with a chain of several '->' ");
            --return null; 
         );
     );

     if numparmsPos===null then 
     (
             bblog.warning (" warning: did not find the position of 'numparms:' in dissasseble string ");
             return null; 
     )
     else
         return  value lst#(numparmsPos+1);
);



-- find out for a method , how many parameters it does accept, only if a single function is installed for that method
-- if multiple functions are installed for the same method, returns null.
-- status: beta.
guessAcceptedParameterNumber( MethodFunction ) := ZZ=>(foo)->
(
   func := apply( methods foo , m-> lookup m);
   if #func==1 then 
   (
        return  (guessAcceptedParameterNumber func#0);
   );
   if #func>1 then 
   (
     bblog.info ("did not expect a method with multiple installed functions. ");
     return null;
   )
   else
   (
     error ("guessAcceptedParameterNumber: no functions installed for that method; something is screwed up ");
   );
);



testGuessAcceptedParameterNumber = ()->
(
   a:=null;   b:=null;   c:=null;

   foo := (a)->(5);
   assert(1==guessAcceptedParameterNumber foo);

   bar := (a,b)->(5);
   assert(2==guessAcceptedParameterNumber bar);

   foobar := method();
   foobar(Boolean,Boolean,String) := ZZ => (a, b, c)->5;
   assert( 3==guessAcceptedParameterNumber foobar );

   -- do not know how to check and what to do for this case:
   -- foo = a->(a,b)->(5);

   foo = a->(a,b);

   assert( 1==guessAcceptedParameterNumber foo );

)


TEST ///
 debug BlackBoxIdeals
 idealBlackBoxesProtect()
 testGuessAcceptedParameterNumber()
///


 

-- Estring = new Type of String - goal was to distinguies string keys from symbol keys without changing string visuaulization (net())
-- but that turned out to be impossible, because e.g. operator '#' cannot be overloaded.


--net (String) := Net =>(str)->(   strl := "\""| str| "\"" ;   return stack separate strl; );

savedEpsRings := new MutableHashTable;

getPropertySymbols := method ();
   getPropertySymbols(String) := List => (propertyName)->
   (
      --
      propertySymbols := {} ;
      --
      try  (  propertySymbols = propertySymbols | { getGlobalSymbol(BlackBoxIdeals.Dictionary, propertyName)} );
      
      -- todo question: should the symbol in the users private dictionary always be created?
      try  (  propertySymbols = propertySymbols | { getGlobalSymbol propertyName} ) else 
       ( 
              propertySymbols =  propertySymbols | { getGlobalSymbol( User#"private dictionary" propertyName); }
        );
      return propertySymbols;
   );



-- package-global symbol for 
geps := getSymbol "eps"; 

getEpsRing = method();

getEpsRing(Ring, ZZ) := Ring => (coeffring, epsDim)->
(
    leps := geps;
 
    if (epsDim<0) then error("expected epsDim>0 ");
    
    epsRng := null;
    if not (savedEpsRings#?(coeffring,epsDim) ) then 
     (
        polRing:=coeffring[leps];
        leps=(gens(polRing))#0;
        savedEpsRings#(coeffring,epsDim) = polRing/leps^(epsDim+1);    
        epsRng = savedEpsRings#(coeffring, epsDim);
        eps = (gens epsRng)#0;
        for symb in getPropertySymbols("eps") do 
        (
           assert(symb=!=null);
           (savedEpsRings#(coeffring,epsDim))#symb  = eps;
        )
    ); 
    epsRng = savedEpsRings#(coeffring, epsDim);
    eps = (gens epsRng)#0;
    return epsRng
)

testEpsRing= ()->
(
    rng := getEpsRing(ZZ,1);
    eps:= (gens rng)#0;
    assert (1_rng+eps_rng + eps_rng*eps_rng == 1_rng+eps_rng );

    rng1 := getEpsRing(ZZ,1);
 
    rng2 := getEpsRing(ZZ,1);

    assert(rng1===rng2 );
);


TEST ///
    debug BlackBoxIdeals
    testEpsRing();
///



-- polynomialLCMDenominator computes the least common multiple of the denominators of the polynomial (rational) cofficients.
-- that means if we have a polynomial = sum { (a_i/b_i)*monomial_i }, a_i and b_i integers, then the function returns LCM( (b_i) )
polynomialLCMDenominator = (polynomial)->
(
    coeffRng := null;
    LCMDenominator := 1;
    --
    summands := { polynomial };
    while (coeffRng=!=ZZ and coeffRng=!=QQ) do
    (
        try ( coeffRng = coefficientRing ring  summands#0 ) then
        (        
             summands = flatten apply( summands, summand-> apply(flatten entries  last coefficients summand, j->sub(j,coeffRng) ) );    
        )
        else
        (  error("expected rationals as coefficient ring!"); 
        );
    );
    if (coeffRng===QQ) then   LCMDenominator =  lcm apply(summands ,j-> denominator j ) ;
    return LCMDenominator;
)



-- clearCoeffDenominators converts an ideal with rational coefficients to an ideal with integer coefficients while preserving the vanishing set.
-- e.g. if sub(IdealWithRationalCoeffs,point)==0, then  sub( clearCoeffDenominators(IdealWithRationalCoeffs),point)==0 and vice versa

clearCoeffDenominators  = method();

clearCoeffDenominators (Ideal)  :=  Ideal =>  (IdealWithRationalCoeffs)->
(
    if (coefficientRing ring IdealWithRationalCoeffs=!=ZZ and coefficientRing ring IdealWithRationalCoeffs=!=QQ) then
    error("expected rationals as coefficient ring!");
    dstrng := ZZ[gens ring IdealWithRationalCoeffs];
    modgens := apply(flatten entries gens IdealWithRationalCoeffs, i->polynomialLCMDenominator(i)*i );
    return sub(ideal modgens,dstrng );
)


doc ///
    Key
        clearCoeffDenominators
        (clearCoeffDenominators, Ideal )   
    Headline
        convert an ideal with rational coefficients to an ideal with integer coefficients
    Usage   
        clearCoeffDenominators(IdealInQQ)
    Inputs  
        IdealInQQ:Ideal
             ideal with rational coefficients
    Outputs
        : Ideal
             ideal with integer coefficients with the same zero set over QQ as the input ideal
    Description
        Example          
        Text
           \break  Example:  convert an ideal with coefficients in QQ to an ideal with   coefficients in ZZ
        Example          
            RQ = QQ[x];
            FQ = {1/3*x+1,1/5*x+2};        
            IFQ = ideal FQ
            IFZ = clearCoeffDenominators(IFQ)
    Caveat
        Conversion implemented only for cases where the ideal coefficient ring is QQ( or ZZ).
///



--        clearCoeffDenominators        

-- test conversion of an ideal with rational coefficients to an ideal with integer coefficients while preserving the vanishing set. 
testClearCoeffDenominators =()->
(
    x := null;  x=symbol x;
    y := null;  y=symbol y;
    RQ := QQ[x,y];
    x = (gens(RQ))#0;
    FQ := { (1/3*x+1) ,  (y+1/2)}; 
    IFQ := ideal FQ;
    IFZ := clearCoeffDenominators(IFQ);  
    FZ := (entries (gens IFZ)_0)#0;

    assert(   FZ == sub(x+3,ring FZ)   ) ; 

    point := matrix {{-3,-1/2}};
    assert( sub(IFQ,point)==0 );
    assert( sub(IFZ,point)==0 );

    point = random(QQ^1,QQ^2);
    assert( sub(IFQ,point)==sub(IFZ,point) );    
)

-- test polynomialLCMDenominator ( computing the least common multiple of the denominators of polynomials with rational coefficients)
-- in case that the rationals(QQ) are not the coefficient ring

testNestedRingCoeffsLCMDenominator =()->
(
    x:=null; x=symbol x;
    y:=null;  y=symbol y;
    z:=null;  z=symbol z;
    RQ := QQ[x,y];
    x = (gens(RQ))#0;

    RQQ := RQ[z];
    polFQ :=  (1/3*x+1)*z; 

    lcmDenom := polynomialLCMDenominator( polFQ );
    assert(lcmDenom==3);   
)

-- test conversion of an ideal with rational coefficients to an ideal with integer coefficients while preserving the vanishing set. 
-- (clearCoeffDenominators)
-- the special case, that the ring the equations belong to is a tensor product of other rings. 
testTensoredClearCoeffDenominators =()->
(
    x:=null; x=symbol x;
    y:=null;  y=symbol y;
    z:=null;  z=symbol z;
    R1Q := QQ[x,y];
    x = (gens(R1Q))#0;

    R2Q := QQ[z];

    RTQ := R1Q**R2Q**QQ;

    x = (gens(RTQ))#0;
    z = (gens(RTQ))#2;
    polFTQ :=  (1/3*x+1)*z; 

    lcmDenom := polynomialLCMDenominator( polFTQ );
    assert(lcmDenom==3);
    IFQ := ideal polFTQ;
    IFZ := clearCoeffDenominators(IFQ);  
    FZ := (entries (gens IFZ)_0)#0;
    assert(   FZ == sub(x*z+3*z,ring FZ)   ) ;       
)

-- deduceNumGenerators:
--   for an (ideal) blackbox, determine the number of   generators (or equations). This is possible in case 
--   the blackbox provides the 'valuesAt'-property.
--
--  todo: this is a generic version. if the blackbox is given by an ideal I, the generators can be determined easily by #(gens ideal I)
--  todo: what should the user do, if deduceNumGenerators fails? 
--
deduceNumGenerators := (blackBox)->
(
     bblog.debug(" enter deduceNumGenerators ") ;

     if not blackBox.hasPointProperty("valuesAt") then error ("cannot determine the number of generators/equations, since valuesAt is missing");
     try ( numVar:=blackBox.numVariables; ) else (  error (" blackBox.numVariables is missing"); );
    
     computed  := false; 
     maxTrials := 100;
     currTrial := 0;
     rng := blackBox.coefficientRing;
     numGenerators := null;
     while numGenerators===null and currTrial<maxTrials do
     (
          try (
 
          bblog.debug(" deduceNumGenerators: enter try block ") ;

              tmppoint := matrix random(rng^1,rng^(blackBox.numVariables) );

              bblog.debug(" deduceNumGenerators:computed tmppoint ") ;

              valuesMatrix := blackBox.valuesAt( tmppoint );
  
               bblog.debug(" valuesMatrix computed ") ;

              --print valuesMatrix;

              assert (numRows valuesMatrix==1);
              numGenerators = numColumns valuesMatrix;

          )         then ( computed=true ) else();
               
           
          currTrial = currTrial+1;
    );
    if not computed then error "error: failed to deduce number of generators/equations; please report the problem to the package developers";
    return numGenerators;
);


testDeduceNumGenerators=()->
(
  blackBoxDummy := new MutableHashTable;
  blackBoxDummy.hasPointProperty = (propertyName)->
  (
     if propertyName==="valuesAt" then return true;
     return false;
  );
  blackBoxDummy.valuesAt= (point)-> ( return matrix {{1,2,3,4,5}} );

  blackBoxDummy.coefficientRing = ZZ;
  blackBoxDummy.numVariables = 5;
  numGenerators := deduceNumGenerators(blackBoxDummy);
  assert( numGenerators==5 );
);


TEST ///
 debug BlackBoxIdeals
 idealBlackBoxesProtect()
 testDeduceNumGenerators()
///


-- dropDegreeInfo:
--   for some matrix operations (which ones?), degree information needs to be dropped,
--    which is done by this method. Used in '.jacobianAt' which in turn is  used in the 'padicLift' package.
--
dropDegreeInfo = method();
dropDegreeInfo (Matrix) := Matrix=> (mat)->
(
   return map( (ring mat)^(numRows mat) ,(ring mat)^(numColumns mat), mat );
);


-- constructs a jacobian at a point supposing that the blackBox implements evaluation at a point.
-- currently expects that blackBox implements 'valuesAt' and knows number of generators/equations.(numGenerators).
-- the second dependency could be removed, as the column number  of the returned 'valuesAt' evaluation ( a row vector) should be the same as number of generators/equations.
-- hmm, this stuff is very hard to debug because of the try clauses in the black box.
--
deduceJacobianAt = ( blackBox, point )->
(
  
    --valuesAt := null; numGenerators := null;
    --valuesAt = global valuesAt;

    --numGenerators = global numGenerators;

    if not  blackBox.hasPointProperty("valuesAt") 
       then error("deduceJacobianAt: to construct jacobian at a point the black box need at least the property 'valuesAt' ");
    -- assert( blackBox#?(global valuesAt) ); -- this should be a 
    try {  blackBox.numGenerators() } else (
       error("deduceJacobianAt: to construct jacobian at a point the black box requested number of generators ('numGenerators'), which was not present.");
    );

    -- assert( blackBox#?(global numGenerators) ); -- this should be a test.

    rngPoint := ring point;
    numVariables := blackBox.numVariables ;      
    if (not ( blackBox.valuesAt( point )==0))  then  error("point does not belong to the ideal ! ");
    
    epsRng := getEpsRing(rngPoint, 1);
    eps := (gens epsRng)#0;


    jacobianMatrixAt := mutableMatrix( rngPoint ,  numVariables, blackBox.numGenerators() );
    for unknownIdx  in 0..(numVariables-1) do
    (
         newpoint := new MutableMatrix from sub( point, epsRng );

         newpoint_(0,unknownIdx) = newpoint_(0,unknownIdx)+eps;
         valueVec := blackBox.valuesAt( matrix newpoint );  
         for equationIdx in 0..numColumns valueVec-1 do
         (
            coordinateValue := last coefficients (valueVec_(0,equationIdx), Monomials=>{1 , eps } );
            if ( not (coordinateValue)_(0,0) ==0) then error("error in jacobianAt. please contact the developers");
            jacobianMatrixAt_(unknownIdx,equationIdx) = sub( (coordinateValue)_(1,0) , rngPoint  )  ;
         );
    );
    return matrix jacobianMatrixAt;
);


testDeduceJacobianAt = ()->
(
  K := ZZ/5;
  x := getSymbol "x";  y := getSymbol "y";  z := getSymbol "z";
  
  R := K[x,y,z];

  x = (gens R)#0;  y = (gens R)#1;  z = (gens R)#2;
 
   -- ideal of a plane and a line that intersect at the origin
   I := ideal (x*z,y*z);

  blackBoxDummy := new MutableHashTable;

  point := matrix{{1_R,1_R,0_R}};

   blackBoxDummy.numGenerators= ()-> return (numColumns (sub (gens I, point ) ) ) ;


  blackBoxDummy.hasPointProperty = (propertyName)->
  (
     if propertyName==="valuesAt" or propertyName==="jacobianAt"  then return true;
     return false;
  );

  blackBoxDummy.valuesAt= (point)-> ( return sub(generators I, point); );

  blackBoxDummy.numVariables = #(gens R);



  computedJac := deduceJacobianAt( blackBoxDummy, point);

  targetJac := dropDegreeInfo sub(jacobian generators I, point); 
  assert( computedJac == targetJac );

  point= matrix{{ 0_R, 0_R, 5_R}};
  computedJac = deduceJacobianAt( blackBoxDummy, point);
  targetJac = dropDegreeInfo sub(jacobian generators I, point); 
  assert( computedJac == targetJac );
);


TEST ///
 debug BlackBoxIdeals
 idealBlackBoxesProtect()
 testDeduceJacobianAt()
///






-- introduce a new type representing a parameter space.

BlackBoxParameterSpace = new Type of  HashTable;



-- jetAtSingleTrial: 
--   tries once to compute a jet , ( see http://en.wikipedia.org/wiki/Jet_%28mathematics%29 for jet definition;) 
--   for the used computation algorithm see the bacherlor thesis at 'http://www.centerfocus.de/theses/js.pdf' .
--
--   preconditions: black box provides evaluation at a point. ('valuesAt') 
-- 
-- returns a hashtable with entries
-- 
--  - "succeeded" a boolean, 
--  - "failedJetLength"  contains the jet length at which the computation failed, otherwise null
--  - "jet"  contains the jet . The jet of the length n has the form
--            j = point + eps*y_1 + . . . + eps^n * y_n 
--           such that F(j) = 0, 
--           where F: E_(n+1)^m -> E_(n+1)^k 
--             with E_(n+1) = K[eps]/( eps^(n+1) ) 
--             whereby K is the coefficient ring (blackBox.coefficientRing), 
--             m is the number of variables (blackBox.numVariables) of the parameter space (same as entries in the point vector)
--             and k is the number of the generators/equation of the (implicitly or explicitly) given ideal. 

--
jetAtSingleTrial = method();

-- here we improve precision by 1 in each step
-- using newtons-Algorithm one could double precision in each step, but
-- for this we also need high precision jacobi-Matrices
-- For black-Box-Jacobi-Matrices we might not have high precision Jacobi-Matrices

--jetAtSingleTrial( BlackBoxParameterSpace, Matrix, ZZ ) := MutableHashTable => ( blackBox,  point, jetLength )  ->
jetAtSingleTrial( HashTable, Matrix, ZZ ) := MutableHashTable => ( blackBox,  point, jetLength )  ->
(

    -- if not (blackBox.isZeroAt(point)) then error(" function is not zero at the point ");
    retVal := null;
    if not (blackBox.isZeroAt(point))  then
    (
       retVal = new HashTable from { "succeeded" => false, "jet" => null, "failedJetLength" =>0 };
       return retVal;
    );

    liftingFailed := false;
    failedJetLength := null;
    jet := point;

    jacobianM2Transposed := transpose blackBox.jacobianAt(point) ;
    
    jtColumns := numColumns jacobianM2Transposed ;
    jtRows    :=  numRows jacobianM2Transposed ;

    jacobianKernel := generators kernel jacobianM2Transposed ; -- syz would also work

    epsRng := getEpsRing( blackBox.coefficientRing,  1 );
    eps = (gens epsRng)#0;

    coeffRng := (blackBox.coefficientRing); -- klammern sind Notwendig!

     
    rnd := random( coeffRng^(numColumns(jacobianKernel)), coeffRng^1 );
    if (numColumns(jacobianKernel)>0) then 
    (   
        while  zero(rnd) do
        (
            rnd = random( coeffRng^(numColumns(jacobianKernel)), coeffRng^1 );
        );
    );
    
    lengthOneLift := sub(point,epsRng) + transpose(sub( jacobianKernel*rnd, epsRng) *eps);

    jet = lengthOneLift;

    epsPrecision := 1;     -- first lift should succeed for a smooth point! 
    lastJetNr    := epsPrecision ; 
    
 

    for  epsPrecision in 2..jetLength do 
    (
         epsRng = getEpsRing( coeffRng, epsPrecision);
        eps = (gens epsRng)#0;
  
        valuesAtJet := blackBox.valuesAt( sub(jet,epsRng) );

        rightHandSide := matrix mutableMatrix( coeffRng, numColumns valuesAtJet ,1 );
        
        -- notwendig:
   
        if not zero(valuesAtJet) then 
            rightHandSide = transpose last coefficients (valuesAtJet, Monomials=>{ eps^epsPrecision });
          -- one could also use contract since eps^epsPrec is the highest possible degree
 
        rightHandSide = sub(rightHandSide,coeffRng);
    
        if not (0==rightHandSide % jacobianM2Transposed ) then (
           failedJetLength = epsPrecision;
           liftingFailed=true; break; 
        );

        x := rightHandSide // jacobianM2Transposed ;
        x = x + jacobianKernel* random(coeffRng^(numColumns(jacobianKernel)), coeffRng^1 );
           x = transpose x;
 
        jet2 := sub (jet, epsRng )-sub( x, epsRng ) *eps^epsPrecision;
        assert ( 0 == blackBox.valuesAt(jet2) );
        jet = jet2;
    );
    -- todo: create a datatype or a hashTable for the return value.

    retVal = new HashTable from { "succeeded" => not liftingFailed, "jet" => jet, "failedJetLength" =>failedJetLength };
    return retVal;
)



-- 'jetAt'
--   tries to compute a jet at a given point (which belongs to the zero set of the black box ) several times (= numTrials). 
-- 
--   for the returned result and input restrictions see 'jetAtSingleTrial'

jetAt = method();

jetAt( BlackBoxParameterSpace, Matrix, ZZ, ZZ) := MutableHashTable => ( blackBox,  point, jetLength, numTrials )  ->
(
    assert( numTrials>=1 );
    bestJet := jetAtSingleTrial ( blackBox,  point, jetLength);
    for i in 2..numTrials do
    (
        if (bestJet#"succeeded") then 
        (
            return bestJet;
        )
        else
        (
            jetTrialResult := jetAtSingleTrial ( blackBox,  point, jetLength);
            if ( bestJet#"failedJetLength"< jetTrialResult#"failedJetLength" ) then
            ( 
                bestJet = jetTrialResult;
            );
        );
    );
    --return  new HashTable from { "succeeded" =>false, "lift" => null };
    return bestJet;
)

-- jetAtWrapper currently deprecated !
-- options makes more sense for trials (=1)
jetAtWrapper = method( Options=>{"jetLength" =>4, "trials"=>5} );

jetAtWrapper( BlackBoxParameterSpace, Matrix )  := MutableHashTable => opts -> ( blackBox,  point )  ->
(
    return jetAt( blackBox,  point, opts#"jetLenth", opts#"trials");
);


isCertainlySingularAt = method();
--isCertainlySingularAt( BlackBoxParameterSpace, Matrix, ZZ, ZZ) := MutableHashTable => ( blackBox,  point, jetLength, numTrials ) ->
isCertainlySingularAt( HashTable, Matrix, ZZ, ZZ) := MutableHashTable => ( blackBox,  point, jetLength, numTrials ) ->
(
    assert( numTrials>=1 );

    for i in 1..numTrials do
    (
        worstJet := jetAtSingleTrial ( blackBox,  point, jetLength);
        if (not worstJet#"succeeded") then 
        (
            return true;
        );
    );
    --return  new HashTable from { "succeeded" =>false, "lift" => null };
    return false; -- should return unknown or something similar!
);

-- 'keysWithoutSymbols': returns hashtable keys without symbol keys.
--
keysWithoutSymbols = method();
keysWithoutSymbols(HashTable) := List => (bb)->
(
   --print "my keys";
   bbh:= new HashTable from bb;
   kl := keys bbh;
   klRes := kl;
   for key in kl do 
   (  
      if (class key)===Symbol then
      (
        strKey := toString key;
        if bbh#?strKey then
              klRes=delete(key,klRes);
      );
   );
  return klRes;
)




-- disable BlackBoxParameterSpace construction using new for an arbitrary parameter:
new BlackBoxParameterSpace from Thing := ( E, thing) -> (
   error "creating blackbox from  type " | toString E | " not implemented ";
);

-- disable BlackBoxParameterSpace construction using new without parameters:
new BlackBoxParameterSpace   := (E) -> (
   error "creating empty blackbox not possible. You have at least to provide the number of variables and their ring";
);


BlackBoxIdeal = new Type of  BlackBoxParameterSpace;


-- disable BlackBoxIdeal construction using new without parameters ( bb = new BlackBoxIdeal ):
new BlackBoxIdeal   := (E) -> (
   error "creating empty blackbox not possible. You have at least to provide the number of variables and their ring";
);

-- disable BlackBoxIdeal construction using new for an arbitrary parameter. ( bb=new BlackBoxIdeal from ... ):
new BlackBoxIdeal from Thing := ( E, thing) -> (

   error "creating blackbox from  type " | toString E | " not implemented ";

);





-- internal method to create a basic black box ( a black box for a parameter space )
--
-- The reason for using an internal (mutable black box) and a public protected blackbox 
-- is to prevent the user from accidental object changing.

-- The object write protection has to be used as  follows: 
--  as a final object for the user always a nonmutable object copy has to be returned;
--  while for the internal 'class' inheritance the non-copied original mutable object (as created by the internal methods) is needed.
--
-- Since the user has access only to the copy and may modify the object e.g. by registering properties,
-- all (potential) mutable black box properties needs to be stored as local variables in the internal methods, 
-- like 'pointProperties' in 'blackBoxParameterSpaceInternal'.
--
-- Access to potential mutable  black box properties by a user can only be modelled  through 'get'-methods, like blackBox.numGenerators()
-- and may never accessed directly, because then other (shallow) copies of the same black box would run out of sync.
-- Because blackBox.numGenerators() is defined in the same context with the local variable 'localNumGenerators', the variable 'localNumGenerators'
-- is visible (and modifiyable) inside 'blackBoxParameterSpaceInternal', but not outside!


blackBoxParameterSpaceInternal = method();

blackBoxParameterSpaceInternal( ZZ, Ring ) := HashTable => ( numVariables, coeffRing ) ->
(
    
   blackBox := new MutableHashTable;   

   -- public: 
   blackBox.coefficientRing = coeffRing;  
   blackBox.numVariables = numVariables;  -- stores the number of the variables in the parameter space.

   -- private: 
   pointProperties := new HashTable  ;     -- stores the point properties
   localNumGenerators := null;                  -- stores the number of the generators/equations. Since this value may be updated,
                                           -- the value cannot be stored in the blackBox directly as a key-value pair, 
                                           -- because, otherwise different black box references referring to the same object
                                           -- could get out of sync. The variable is accessed by a getter(numVariables())

    singularTestOptions := new MutableHashTable;
    singularTestOptions.precision = 10;
    singularTestOptions.numTrials = 2;

   
     -- checks the consistency of the point with the black box 
     -- 2. the point should be given as a column matrix having same number of columns as blackBox.numVariables
     -- 1. if blackBox.coefficientRing is ZZ, then every ring for point entries is allowed (a) why, , (b) is this always correct?
     --    if blackBox.coefficientRing is NOT ZZ, then either the ring of point has to coincide with the 'blackBox.coefficientRing', or
     --                                             the coefficientRing( ring point) has to coincide with the 'blackBox.coefficientRing'. 
     -- 
     checkInputPoint := (point)->
     (
        errorMsg := "  ideal is defined over "| toString  blackBox.coefficientRing | "and not over " | toString ring point |"!" ;
        if blackBox.coefficientRing =!= ZZ  
           and ring point =!= blackBox.coefficientRing
           then 
       (
           try ( if coefficientRing ring point=!= blackBox.coefficientRing then error errorMsg )  then () else ( error (errorMsg) );
           );
       if (numColumns point=!=blackBox.numVariables) then error (" illegal point : expected " | toString blackBox.numVariables | " coordinates");
     );


     -- pointProperty:
     --    return a point property (by symbol or by name) stored in 'pointProperties'. 
     --    If a property does not exist, throws an error.
     --
     blackBox.pointProperty = method();

     blackBox.pointProperty (String) := Function =>( prop ) ->
     (
        if not pointProperties#?prop then        
            error (" blackbox does not have requested point property " | toString prop  );

        return pointProperties#prop;
     );     
     
     blackBox.pointProperty (Symbol) := Function =>( prop ) ->
     (
        return blackBox.pointProperty(toString prop );
     );


   --  hasPointProperty:
   --    checks, if the blackbox has a specified property (by name or by symbol)
   --
   blackBox.hasPointProperty = method()  ;
   blackBox.hasPointProperty (String) := Boolean =>(propertyName)->
   (
       return pointProperties#?propertyName;
   );

   blackBox.hasPointProperty ( Symbol ) := Boolean => (propertySymbol)->
   (
       return blackBox.hasPointProperty(toString propertySymbol)
   );

    


   -- 'knownPointProperties':   
   --     returns a list of all registered point properties (as Strings)
   --
   blackBox.knownPointProperties = ()->
   (   
      return sort unique apply (keys pointProperties, key->toString key);
   );


   -- 'knownPointPropertiesAsSymbols' : 
   --   returns a list of all registered point properties (as Symbols)
   --



   getPropertySymbol := method ();
   getPropertySymbol(String) := Symbol=> (propertyName)->
   (
     
      propertySymbol := null;
      try  (  propertySymbol = getGlobalSymbol propertyName; ) else ( 
              propertySymbol = getGlobalSymbol( User#"private dictionary", propertyName); 
       );
      -- todo first question : is the behaviour above same as for 'global 'symbol'? 
      -- second question: should 
      return propertySymbol;
   );

   blackBox.knownPointPropertiesAsSymbols = ()->
   (   
      return apply( blackBox.knownPointProperties(), propertyName-> getPropertySymbol(propertyName) );
   );


   



   -- 'setPointProperty': internal method to set a point property. 
   --                     Is called by 'outerSetPointProperty', 'setIsZeroAt', 'setValuesAt', 'setJacobianAt'
   -- 
   -- the method works as follows: 
   -- 1. the current (internal) variable 'pointProperties' is copied and transformed to a mutable HashTable
   -- 3. the propertyMethod (see 2) is added to the internal variable  'pointProperties'
   -- 4. the 'pointProperties' are changed to immutable.
   -- 5. the (internal) blackBox HashTable is extended by methods which accept one parameter and
   --     call the corresponding propertyMethod in 'pointProperty' This level of indirection is done,
   --     to keep access to the correct pointProperties even if the internal BlackBox object is replaced or updated: 
   --     all (updated) blackboxes will refer to the same 'pointProperties' variable. 
   --  

   --  after a call of 'getUpdatedBlackBox', the property is accessible by its symbol name and eventually by the symbol, if there is no symbol clash
   -- 
   -- Remark. the first parameter is a symbol and not a name, because it is imaginable, that a user / a package author 
   --   could want to pass the concrete symbol he wants.  
   --
   setPointProperty := method();
   setPointProperty ( Symbol, Function ) := Thing => ( propertySymbol, propertyMethod )->
   (

     bblog.debug(" called setPointProperty ") ;

      -- propertySymbol := getGlobalSymbol propertyName;
      propertyName := toString propertySymbol;

      assert(propertyName=!=null); 
      
      packageSymbol := null;
      if BlackBoxIdeals.Dictionary#?propertyName then 
        packageSymbol  = BlackBoxIdeals.Dictionary#propertyName;
    
      -- step 1
      pointProperties = new MutableHashTable from pointProperties;

     bblog.debug(" set point Property " | propertyName) ;

 
     -- step 2,3
          -- remember: cannot return in a try statement!
          pointProperties#propertyName = ( point )-> ( 
          checkInputPoint(point);
            result := (propertyMethod)(  point ); 
          return result;
     );

     -- inconsistent and unnecessary ?    
     -- pointProperties#propertySymbol = pointProperties#propertyName;
  
     if packageSymbol=!=null then 
       pointProperties#packageSymbol = pointProperties#propertyName;

  
  
     -- step 5
     if mutable blackBox then 
     (
         blackBox#propertySymbol        = (point)->( (blackBox.pointProperty(propertyName))(point) );
         blackBox#propertyName          = (point)->( (blackBox.pointProperty(propertyName))(point) );

         --if packageSymbol=!=null then 
         --  blackBox#packageSymbol        = (point)->( (blackBox.pointProperty(propertyName))(point) );

         for symb in getPropertySymbols(propertyName) do 
         (
             assert(symb=!=null);
             blackBox#symb  = (point)->( (blackBox.pointProperty(propertyName))(point) );
         )
     );

      -- step 4
      pointProperties = new HashTable from pointProperties;
  
   );

   setPointProperty ( String, Function) := Thing => ( propertyName, propertyMethod )->
   (
       propertySymbol := getPropertySymbol(propertyName);
       setPointProperty( propertySymbol, propertyMethod );
   );


   --  valuesAtWrapper: post check the result of the 'valuesAt' method:
   --                  the result is expected to be an 1-row matrix. 
   --
   valuesAtWrapper := ( pValuesAt, point) ->
   (
     result := pValuesAt(  point); 
       try  ( assert(numRows result==1); ) else { error ( "'valuesAt' did not return an one-row matrix "); };
      return result;
   );

   -- setIsZeroAt
   --   sets the check, if a point belongs to the object (isZeroAt(point)==true ) or not.
   --
   --   called by   'outerSetPointProperty'<-{'registerPointProperty', 'updatePointProperty'}, 'setValuesAt'
   -- 
   setIsZeroAt := (pIsZeroAt) ->
   (   
       -- parameter is called differently to the symbol 'isZeroAt', otherwise it seems we could get the wrong value...
       setPointProperty("isZeroAt" , pIsZeroAt );
   );


   -- setJacobianAt:
   --    set a method to compute the jacobian at a point.
   --
   --   called by 'outerSetPointProperty'<-{'registerPointProperty', 'updatePointProperty'}, 'setValuesAt'
   --  
   --   triggers updates for , 'rankJacobianAt'.
   --  
   setJacobianAt := (pJacobianAt) ->
   ( 
      bblog.info( "setJacobianAt: updates also (   rankJacobianAt)" );      

      localJacobianAt := ( point)->
      (
         return dropDegreeInfo( pJacobianAt(point) );
      );

      setPointProperty( "jacobianAt" , localJacobianAt );

      localRankJacobianAt := (  point)->
      (
         return  rank blackBox.jacobianAt(point) ;
      );

      localIsCertainlySingularAtWrapper  := ( point )  ->
        (
               return isCertainlySingularAt( blackBox,  point, singularTestOptions.precision,  singularTestOptions.numTrials );
        );

      setPointProperty("isCertainlySingularAt" , localIsCertainlySingularAtWrapper );

      localIsProbablySmoothAtWrapper  := ( point )  ->
        (
               return not isCertainlySingularAt( blackBox,  point, singularTestOptions.precision,  singularTestOptions.numTrials );
        );

      setPointProperty("isProbablySmoothAt" , localIsProbablySmoothAtWrapper );


      setPointProperty( "rankJacobianAt" , localRankJacobianAt );
   );


   -- setValuesAt:
   --    set a method to compute the values of the generators/equations at a given point.
   --
   --   called by 'outerSetPointProperty'<-{'registerPointProperty', 'updatePointProperty'}, 'setValuesAt'
   --  
   --   triggers updates for 'isZeroAt', 'numGenerators', jacobianAt', 'rankJacobianAt'.
   --  
   setValuesAt := (pValuesAt) ->
   (      
       bblog.info( "setValuesAt: updates (isZeroAt, numGenerators, jacobianAt, rankJacobianAt)" );      

       localValuesAt := (point)->return valuesAtWrapper(pValuesAt, point ) ;

       -- when using valuesAt instead of localValuesAt we get the wrong  (symbol valuesAt) (local valuesAt)

       setPointProperty( "valuesAt"  ,  localValuesAt );
       
        localNumGenerators =  deduceNumGenerators(blackBox)  ; --depends on valuesAt.

       blackBox.numGenerators = ()->(return localNumGenerators);


       bblog.info( "updated blackBox.numGenerators to " | toString blackBox.numGenerators() );   
     


       setIsZeroAt(
          (point)->( return blackBox.valuesAt(point)==0 ;) 
        );

       ----- jacobian at:
       localMethod :=  (point)->deduceJacobianAt( blackBox, point );
       setJacobianAt ( localMethod );  
   );
 

   -- outerSetPointProperty:
   --    this method is the common between 'updatePointProperty' and 'registerPointProperty' and was therefore outsorced
   --   setting properties isZeroAt, valuesAt, jacobianAt are handled especially.

   --    since it is allowed for the provided propertyMethod call to accept one parameter (point ) 
   --      or two parameters (blackBox, point) the two possible calls are wrapped in a function, which accepts a single parameter (point)

   outerSetPointProperty := ( propertySymbol, propertyMethod)->
   (
     propertyName := toString  propertySymbol;


      acceptedNumParameters := guessAcceptedParameterNumber propertyMethod;

      if not (acceptedNumParameters==2 or  acceptedNumParameters==1 ) then 
          error (" provided method " | propertyName | " expected to accept 1 or 2 parameters:  ( blackbox, point ),  or (point) , but the passed one seems to accept " | toString acceptedNumParameters);



      --if acceptedNumParameters=!=2 then 
      --    error (" provided method " | propertyName | " expected to accept 2 parameters ( blackbox, point ),  but the passed one seems to accept " | toString acceptedNumParameters);
      
       -- now wrap the provided method if neccesary in a way that it accepts only a point: 
 
       localPropertyMethod := propertyMethod;

        if acceptedNumParameters==2 then 
          localPropertyMethod = ( point )-> ( 
             return propertyMethod( blackBox,   point ); 
        );
  
    

      
         if propertyName==="isZeroAt" then 
              return setIsZeroAt(localPropertyMethod); --probably not necessary

         if propertyName==="valuesAt" then 
             return setValuesAt(localPropertyMethod);  -- triggers initialization of 'isZeroAt' , 'numGenerators' and 'jacobianAt'
     
         if propertyName==="jacobianAt" then 
             return setJacobianAt(localPropertyMethod); -- triggers initialization of  'rankJacobianAt'

           setPointProperty( propertySymbol, localPropertyMethod );
   );

   -- todo : test ; three scenarios should work: 
   --          a user registers a point property
   --          a point property is registered in this package
   --          a point property is registered in a different package
   -- 
  


   -- todo: test if updating 'valuesAt' will  trigger updating isZeroAt and jacobianAt 
   -- 
  
   --  updatePointProperty:
   --    if a property is already set, updates it.
   --
   blackBox.updatePointProperty = method();

   blackBox.updatePointProperty(String, Function) := (BlackBoxParameterSpace) => ( propertyName, propertyMethod )->
   (
      if  (  pointProperties#?propertyName ) then
      (  
           propertySymbol := getPropertySymbol(propertyName);
           outerSetPointProperty( propertySymbol, propertyMethod );
      ) 
      else ( error ("point property "| toString propertyName | " does not exist.") );
      return blackBox.getUpdatedBlackBox(blackBox);
   );



   -- registerPointProperty: a method to register a point property, while providing a propertySymbol,  
   --  expecting that after registering (and getUpdatedBlackBox() ) the property will be accessible via  blackBox#propertySymbol .
   --  usually providing the corresponding symbol is not necessary, but it could be, since each package has its own symbol scope.
   -- todo: this one should not be a public one

   blackBox.registerPointProperty = method();


   blackBox.registerPointProperty(String, Symbol, Function) := BlackBoxParameterSpace => ( propertyName, propertySymbol, propertyMethod )->
   (

      assert( (toString propertySymbol)==propertyName);

      if  ( not  pointProperties#?propertyName 
      and not  blackBox#?propertySymbol          and   not  blackBox#?propertyName  ) then 
      (
           outerSetPointProperty( propertySymbol, propertyMethod );
      )
      else error(" key  "| propertyName |"  exists already. If it is a point property, please use 'updatePointProperty' for updating.");
      return blackBox.getUpdatedBlackBox(blackBox);
   );


   blackBox.registerPointProperty(String, Function) := Thing => ( propertyName, propertyMethod )->
   (
      propertySymbol :=  getPropertySymbol(propertyName);
      return blackBox.registerPointProperty(  propertyName, propertySymbol, propertyMethod )
      --return blackBox.getUpdatedBlackBox(blackBox);
   );

   blackBox.rpp =  blackBox.registerPointProperty;

  blackBox.upp =  blackBox.updatePointProperty;

   -- return a list of known methods. Manually updated.
   blackBox.knownMethods = ()->
   (   
      methods:= { getGlobalSymbol( BlackBoxIdeals.Dictionary, "knownMethods" ) ,
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "knownAttributes" ) ,
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "knownPointProperties" ),
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "knownPointPropertiesAsSymbols" ),
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "hasPointProperty" ),
                  --getGlobalSymbol( BlackBoxIdeals.Dictionary, "pointProperty" ),
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "registerPointProperty" ),
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "setSingularityTestOptions" ),
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "rpp" ), 
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "upp" ), 
                  getGlobalSymbol( BlackBoxIdeals.Dictionary, "updatePointProperty" )
                  ----getGlobalSymbol( BlackBoxIdeals.Dictionary, "getUpdatedBlackBox" ),
                  --getGlobalSymbol( BlackBoxIdeals.Dictionary, "unknownIsValid" )
            };
    --  methods:= {   knownMethods,
    --                knownAttributes,
    --                knownPointProperties,
    --                knownPointPropertiesAsSymbols,
    --                hasPointProperty,
    --                pointProperty,
    --                registerPointProperty, 
    --                updatePointProperty,
    --                getUpdatedBlackBox,
    --                unknownIsValid
    --              };

        sortedMethods := sort apply(methods, i-> ( toString i, i ));
        sortedMethods = apply(sortedMethods, i->(i_1));
        return sortedMethods;  
        --return methods;
   );
  
  
  blackBox.knownProperties = ()->
   (  
      properties := {};
      if ( blackBox#?"valuesAt" ) then
      (
          properties = properties | { getGlobalSymbol( BlackBoxIdeals.Dictionary, "numGenerators" ) };
      );
      return properties;
  );

   -- 'knownAttributes' returns a list of known Attributes. 
   -- (computed as 'keys blackBox which are not Functions' \ {knownPointPropertiesAsSymbols() }
   blackBox.knownAttributes = ()->
   (
         all :=  keysWithoutSymbols blackBox;
         all = select (all, (foo)->(not instance(blackBox#foo,Function)));   
         kM := kP := kPP := kPS := {};
         -- kM =  blackBox.knownMethods();
         -- kP =  blackBox.knownProperties();
         -- kPP =  blackBox.knownPointProperties();
         kPS  =  blackBox.knownPointPropertiesAsSymbols();
         toRemove := kM |  kP | kPP |kPS;
         for symb in toRemove do
           all = delete(symb,all);

         sortedAttributes := sort apply(all, i-> ( toString i, i ));
         sortedAttributes = apply(sortedAttributes, i->(i_1));
         return sortedAttributes;  
     );


    
    
    
    blackBox.setSingularityTestOptions = (prec, numTrials)->
    (
       assert( prec >= 0 and numTrials>0 );
       singularTestOptions.precision = prec;
       singularTestOptions.numTrials = numTrials;
       if blackBox.hasPointProperty("valuesAt") then
       (
           localIsCertainlySingularWrapper  := (  point )  ->
           (
               return isCertainlySingularAt( blackBox,  point, singularTestOptions.precision,  singularTestOptions.numTrials );
           );

           if blackBox.hasPointProperty("isCertainlySingularAt") then 
           (
               blackBox.updatePointProperty("isCertainlySingularAt", localIsCertainlySingularWrapper)
           )
           else
           (
               blackBox.registerPointProperty("isCertainlySingularAt", localIsCertainlySingularWrapper);
           );
           localIsProbablySmoothWrapper  := (  point )  ->
           (
               return isProbablySmoothAt( blackBox,  point, singularTestOptions.precision,  singularTestOptions.numTrials );
           );

           if blackBox.hasPointProperty("isProbablySmoothAt") then 
           (
               blackBox.updatePointProperty("isProbablySmoothAt", localIsProbablySmoothWrapper)
           )
           else
           (
               blackBox.registerPointProperty("isProbablySmoothAt", localIsProbablySmoothWrapper);
           );


      );

    );


   

   blackBox.type = BlackBoxParameterSpace;

   -- a user should not call this method...
   blackBox.getUpdatedBlackBox  = (bb1)->
   (
       
       bb := new MutableHashTable from  bb1;
       for key in keys blackBox do
       (
          if not bb#?key then
          (
              bb#key = blackBox#key;
          )
       );
       blackBox = bb;

       for  property in blackBox.knownPointProperties() do
       (
          propkeys := getPropertySymbols( property ) | {toString property};

          for key in propkeys  do
          (
              --if not bb#?property then 
              if not blackBox#?key then 
              (
                  blackBox.pointProperty(property);
                  blackBox#key = (point)->( (blackBox.pointProperty(toString property))(point) );
              )
              else 
              (      
              );     
         );    

       );

       --bb = new HashTable from bb;
       -- todo: here we have an issue; how to find out the class of the outer blackbox? 
      
       --bb = newClass( BlackBoxParameterSpace, bb );
       bb = newClass(  blackBox.type, bb );
       return  bb;
   );


  
   return blackBox;
)

--            Provide via registerPointProperty() added new black box properties
--          via point('.') operator (syntactic sugar)
rebuildBlackBox = method() ;

rebuildBlackBox(HashTable) := HashTable => ( bb )->
(
       return  bb.getUpdatedBlackBox(bb);
);



-- this function may be used to create a derived object, which inherits properties of an black box ideal
-- since the blackbox object is not copied 
--
blackBoxParameterSpaceInternal(Ring) := HashTable => ( pRing ) ->
(

   blackBox := blackBoxParameterSpaceInternal(#(gens pRing), coefficientRing pRing);

   blackBox.ring = pRing;
   blackBox.unknowns = gens blackBox.ring;
   assert( blackBox.numVariables == #blackBox.unknowns );

   blackBox.unknownIsValid = (unknown)->
   (
        if not ( blackBox.ring === ring unknown) then 
        ( 
            bblog.error( "the unknown is not element of the equations ideal ring" );
            return false;
        );
        return true;
   );


   return blackBox;
)



blackBoxIdeal = method();

 
blackBoxParameterSpace = method();


-- this function is final, that means nobody should use this method for creating a derived object
new BlackBoxParameterSpace from Ring := (E, pRing )->
(
    blackBox := blackBoxParameterSpaceInternal( pRing);
   
    blackBox.type = BlackBoxParameterSpace;
    blackBox = blackBox.getUpdatedBlackBox(blackBox);
    return new HashTable from blackBox;
)


-- this function is final, that means nobody should use this method for creating a derived object
blackBoxParameterSpace(Ring) := HashTable => ( pRing ) ->
( 
    return new BlackBoxParameterSpace from pRing;
);


--blackBoxParameterSpaceInternal(ZZ,Ring) := HashTable => ( numVariables, coeffRing ) ->
--(
--    assert ( numVariables>0 );
--    a := null;
--    a = symbol a;
--    rng := coeffRing[a_1..a_numVariables];
--    return blackBoxParameterSpaceInternal( rng );
--)





-- this function is final, that means nobody should use this method for creating a derived object

blackBoxParameterSpace(ZZ, Ring) := BlackBoxParameterSpace => ( numVariables, coeffRing )  ->
(
    blackBox := blackBoxParameterSpaceInternal( numVariables, coeffRing );
    blackBox.type = BlackBoxParameterSpace;
    blackBox = blackBox.getUpdatedBlackBox(blackBox);
    bb := newClass( BlackBoxParameterSpace, blackBox );
    return bb;
)




-- todo: how to check, if 'ring equationsIdeal' is not a quotient ring?

-- blackBoxIdealInternal
-- this function may be used to create a derived object, which inherits properties of an black box ideal
-- since the blackbox object is not copied 
--
blackBoxIdealInternal := ( equationsIdeal)->
(   
     blackBox :=  blackBoxParameterSpaceInternal( ring equationsIdeal );

    blackBox.type = BlackBoxIdeal;
   
     -- maybe blackBox.addProperty( ideal, equationsIdeal)
     blackBox.ideal = equationsIdeal;      

     -- maybe blackBox.addProperty( equations, gens equationsIdeal )
     blackBox.equations =  gens  equationsIdeal; 


    -- registering ValuesAt generates 'isZeroAt' and 'jacobianAt', too !  
    blackBox.registerPointProperty("valuesAt",( bb, point)->  (    return  gens sub( equationsIdeal , point);   ));


    -- maybe blackBox.addProperty( jacobian, jacobian gens  equationsIdeal )
    blackBox.jacobian = jacobian gens  equationsIdeal;


    -- needs updatePointProperty, because "jacobianAt" is present. 
    -- Todo: introduce force mode option? (overwriting registered property?)
    -- todo: introduce for all methods a precomposition with checkInputPoint
    blackBox.updatePointProperty( "jacobianAt",

       ( bb, point )->
       (   
          jacobianM2MatrixAt := sub( blackBox.jacobian , point);
          return jacobianM2MatrixAt;
       )
    );   
    

    return new HashTable from blackBox; 
)

-- this function is final, that means nobody should use this method for creating a derived object
new BlackBoxIdeal from Ideal := (E, equationsIdeal)->
(
   blackBox := blackBoxIdealInternal(equationsIdeal);
    blackBox = blackBox.getUpdatedBlackBox(blackBox);
    return new HashTable from blackBox;
)

-- this function is final, that means nobody should use this method for creating a derived object
blackBoxIdeal (Ideal) := BlackBoxIdeal =>(equationsIdeal)->
(
   return new BlackBoxIdeal from equationsIdeal;
)




testBlackBoxIdeal=()->
(
   x  := null;
   x  = symbol x;
   rng := ZZ/7[x];
   coeffRng := coefficientRing rng;
   x = (gens rng)#0;

   RP := ZZ/7[x];
   IFP := ideal { 3*x^2+1, 5*x-1 };        
   IFPBlackBox := blackBoxIdeal( IFP );
   point := matrix {{3}};
   rng13 := ZZ/13;
   assert( IFPBlackBox.unknowns=={x} );
   assert( IFPBlackBox.equations==gens IFP);
   assert( IFPBlackBox.jacobian== jacobian IFP);
   
   jac := null;
   point = matrix {{3}};
   try (   jac= IFPBlackBox.jacobianAt(point); ) then 
   (
     error("testblackBoxIdeal: jacobianAt should fail due the coefficient ring of the point matrix does not match the ideal coefficient ring and also the ideal coefficients are not integers ");
   )  else ();
   point = sub( point, coeffRng ) ;
   try {    jac= IFPBlackBox.jacobianAt(point); } else
   (
        error("testblackBoxIdeal: jacobianAt should succeed  due the coefficient ring of the point matrix matches the ideal coefficient ring ");
   );  
    IFPBlackBox.ring;
    IFPBlackBox.valuesAt(point) ;
    assert(   IFPBlackBox.isZeroAt( point ) );
   assert( IFPBlackBox.jacobianAt(point)==sub( jacobian IFP,point) );
   assert( IFPBlackBox.valuesAt(point)== gens sub(  IFP, point ) );
)




-- blackBoxIdealFromEvaluationInternal
--   this function may be used to create a derived object, which inherits properties of an black box ideal
--   since the blackbox object is not copied (and not write-protected)
--
blackBoxIdealFromEvaluationInternal := method();

blackBoxIdealFromEvaluationInternal(ZZ, Ring, Function) := HashTable => ( numVariables, coeffRing, pValuesAt )  ->
(
    blackBox := blackBoxParameterSpaceInternal( numVariables, coeffRing );

    blackBox.type = BlackBoxIdeal;

    blackBox.registerPointProperty ("valuesAt", (bb,point)->pValuesAt(point) ); --sets isZeroAt, jacobianAt, rankJacobianAt and numGenerators

    check := ()->
    (
         numVariables :=  blackBox.numVariables;

         point := matrix { apply(numVariables, i-> 0_(blackBox.coefficientRing) ) };
         blackBox.valuesAt( point );
         blackBox.isZeroAt( point );
    );

    check();  
    return blackBox ;
)


--  creates a BlackBoxIdeal from a given evaluation method ('valuesAt') which takes a point (a row matrix)
--     parameters: numVariables in the parameter space,
--                 coefficientRing , the 
--                 valuesAt: a method for evaluating the object at one parameter point.
--


-- this function is final, that means nobody should use this method for creating a derived object
blackBoxIdealFromEvaluation = method();
blackBoxIdealFromEvaluation(ZZ, Ring, Function) := HashTable => ( numVariables, coeffRing, valuesAt )  ->
(
 
    blackBox := blackBoxIdealFromEvaluationInternal( numVariables, coeffRing, valuesAt ) ;
    blackBox.type = BlackBoxIdeal;
    blackBox = blackBox.getUpdatedBlackBox(blackBox);
    blackBox = newClass( BlackBoxIdeal, blackBox ); 
    return blackBox ;
)


blackBoxIdealFromEvaluation( Ring, Function ) := HashTable => ( pRing, pValuesAt ) ->
(

   blackBox := blackBoxParameterSpaceInternal(pRing );
   blackBox.type = BlackBoxIdeal;
   blackBox.registerPointProperty ("valuesAt", (bb,point)->pValuesAt(point) ); --sets isZeroAt, jacobianAt, rankJacobianAt and numGenerators

    check := ()->
    (
         numVariables :=  blackBox.numVariables;

         point := matrix { apply(numVariables, i-> 0_(blackBox.coefficientRing) ) };
         blackBox.valuesAt( point );
         blackBox.isZeroAt( point );
    );

   check(); 
   blackBox = blackBox.getUpdatedBlackBox(blackBox);
   blackBox = newClass( BlackBoxIdeal, blackBox ); 
   return blackBox;
)




testBlackBoxIdealFromEvaluation = ()->
(
    x  := null;
    x  = symbol x;
    rng := ZZ/7[x];
    coeffRng := coefficientRing rng;
    x = (gens rng)#0;


    RP := ZZ/7[x];
    IFP := ideal ( 3*x^2+1, 5*x-1 );      

    rank source gens IFP;
  
   evaluation := blackBoxIdeal( IFP );
    --evaluation := basicBlackBox();
  
    evalBlackBox := blackBoxIdealFromEvaluation ( # (gens evaluation.ring), coefficientRing evaluation.ring, evaluation.valuesAt );

    point := matrix {{3_(ZZ/7)}} ;
    assert( evaluation.isZeroAt( point ) );
    assert( evaluation.unknowns=={x} );
    assert( evaluation.valuesAt( point ) == evalBlackBox.valuesAt( point ) );
    assert( evaluation.jacobianAt( point ) == evalBlackBox.jacobianAt( point ) );

    assert( evaluation.numVariables ==evalBlackBox.numVariables );

    outerPoint := matrix {{2_(ZZ/7)}} ;

    assert( evaluation.valuesAt( outerPoint ) == evalBlackBox.valuesAt( outerPoint ) );

    assert( not evaluation.isZeroAt( outerPoint ) );

    assert( evaluation.coefficientRing === coeffRng);

    assert( evaluation.unknownIsValid ( ( gens ring IFP )#0 ));

    y  := null;    y  = symbol y;
    rngy := ZZ/7[y];
    y = (gens rng)#0;

    assert( not evaluation.unknownIsValid (  y ) );

);
    

doc ///
   Key
        BlackBoxIdeal
        (NewFromMethod, BlackBoxIdeal, Ideal)
        (NewFromMethod, BlackBoxIdeal, Thing)
   Headline
        an unified interface to a black box ideal
   Usage   
        new BlackBoxIdeal from pIdeal
        blackBoxIdealFromEvaluation( rng, evaluationMethod )
        blackBoxIdealFromEvaluation( variableNumber, coeffRing, evaluationMethod)
   Inputs  
        pIdeal:Ideal
             with integer coefficients
   Outputs
        : BlackBoxIdeal
   Description
         Text
            The {\tt  BlackBoxIdeal } objects implements the interface of @TO BlackBoxParameterSpace @ \break 
            and in addition following methods and attributes:

            \,\, \bullet \,{\tt numGenerators() }: number of  generators/equations \break \break             

            Point properties:\break
            \,\, \bullet \, @TO "isZeroAt" @ \break
            \,\, \bullet \, @TO "valuesAt" @ \break
            \,\, \bullet \, @TO "jacobianAt" @ \break
            \,\, \bullet \, @TO "rankJacobianAt" @ \break
                        
            If the black boxes was generated from an explicit ideal, the
            following properties are also defined: \break 
            \,\, \bullet \,{\tt unknowns}: a list of the parameter variables . \break
            \,\, \bullet \,{\tt ideal}: the original ideal  \break
            \,\, \bullet \,{\tt equations}: generators/equations of the original ideal/polynomial system \break
            
        Text
            \break  For an example see @TO blackBoxIdealFromEvaluation@, @TO blackBoxIdeal @.
///

doc ///
   Key
        blackBoxIdealFromEvaluation        
   Headline
        create a  BlackBoxIdeal  describing an ideal  with integer coefficient ring
   Usage   
        blackBoxIdealFromEvaluation( rng, evaluationMethod )
        blackBoxIdealFromEvaluation( variableNumber, coeffRing, evaluationMethod)
   Inputs  
        rng:Ring
             with integer coefficients
        variableNumber:ZZ
             number of parameter variables
        coeffRing:Ring
             ring to whicht the parameter variables belong to.
        evaluationMethod:Function
             function which accepts a point given as a row matrix and returns the evaluation at this point as a row matrix.
   Outputs
        : BlackBoxIdeal
   Description   
        Text    
           Creates a blackbox describing an implicitly given ideal from an evaluation method \break
           \break  Example:  create an {\tt BlackBoxIdeal } object from an evaluation:
        Example          
            RQ := QQ[x];
            IFQ := ideal { 1/3*x^2-100/3, 1/5*x+2 };        
            evaluation := (point)-> (return sub( gens IFQ, point) );
            bbI := blackBoxIdealFromEvaluation( RQ, evaluation );
        Text
            \break Now access some ideal propeties via the black box interface:
        Example          
            -- keys bbI
            bbI.knownAttributes()
            bbI.unknowns
            point := matrix { {1_QQ } };
            bbI.knownPointProperties()
            bbI.valuesAt(point)

            point = matrix { {-10_QQ } };
            bbI.jacobianAt(point)            
   Caveat
        may have problems if the coefficient ring is something obscure? Rationals and integers as coefficient ring should be ok.
///

--        new BlackBoxIdeal from Ideal
doc ///
   Key
        blackBoxIdeal
   Headline
        create a {\tt BlackBoxIdeal } 
   Usage   
        blackBoxIdeal(anIdeal)
   Inputs  
        anIdeal:Ideal
   Outputs
        : BlackBoxIdeal
   Description      
        Text
           Creates a blackbox describing an ideal \break \break
           See also @TO blackBoxIdealFromEvaluation@  \break
           \break  Example:  create an {\tt BlackBoxIdeal } object from an ideal:
        Example          
            RQ := QQ[x];
            IFQ := ideal { 1/3*x^2-100/3, 1/5*x+2 };        
            IFZ := clearCoeffDenominators(IFQ)
            bbI := blackBoxIdeal(IFZ);
        Text
            \break Now access some ideal propeties via the black box interface:
        Example          
            -- keys bbI
            bbI.knownAttributes()
            bbI.ideal
            bbI.unknowns
            bbI.equations
            bbI.jacobian
            point := matrix { {1_QQ} };
            bbI.knownPointProperties()
            bbI.valuesAt(point)
            point = matrix { {-10_QQ} };
            bbI.valuesAt(point)
            bbI.jacobianAt(point)            
   Caveat
        does not check if the ideal ring is a quotient ring (not supported?)
///


TEST ///
debug BlackBoxIdeals
idealBlackBoxesProtect()
testClearCoeffDenominators()
///



TEST ///
debug BlackBoxIdeals
idealBlackBoxesProtect()
testNestedRingCoeffsLCMDenominator()
///
         

TEST ///
debug BlackBoxIdeals
idealBlackBoxesProtect()
testTensoredClearCoeffDenominators()
///

TEST ///
debug BlackBoxIdeals
idealBlackBoxesProtect()
testBlackBoxIdeal()
///

TEST ///
debug BlackBoxIdeals
idealBlackBoxesProtect()
testBlackBoxIdealFromEvaluation()
///


TEST ///


debug BlackBoxIdeals
idealBlackBoxesProtect()

R = ZZ[x_0..x_3]

M = matrix{
     {x_0,x_1,0},
     {x_1,x_2,x_3}
     }

I = minors(2,M)

smoothPoint = matrix{{1,0,0,0}}
singularPoint = matrix{{0,0,1,0}}
offPoint = matrix{{1,11,0,0}}

S = ZZ[s,t]

line = matrix{{0,0,s,t}}

B = blackBoxIdeal I

assert (sub(jacobian I,smoothPoint) == B.jacobianAt smoothPoint)
assert (rank(B.jacobianAt smoothPoint) == codim(I,Generic=>true))
assert (rank(B.jacobianAt singularPoint) < codim(I,Generic=>true))
assert( B.isZeroAt(smoothPoint) );
assert( B.isZeroAt(singularPoint));
assert(B.isZeroAt(line));

B.valuesAt(offPoint)
B.isZeroAt(offPoint)

assert(not B.isZeroAt(offPoint))

prime = 11;
K = ZZ/prime;
assert (B.isZeroAt(sub(smoothPoint,K)));
assert(B.isZeroAt(sub(offPoint,K)));

apply(100,i->(
      r = random(K^1,K^4);
      assert ((B.isZeroAt r) == ((B.valuesAt r) == 0));
      ));

      
evalLinePlusConic = point -> (
     ePoint = flatten entries point;
     M = matrix{
      {ePoint#0,ePoint#1,0},
      {ePoint#1,ePoint#2,ePoint#3}
      };
      matrix{{det M_{0,1},det M_{0,2},det M_{1,2}}}
     )

B2 = blackBoxIdealFromEvaluation( 4, ZZ, evalLinePlusConic)

apply(100,i->(
      r = random(K^1,K^4);
      assert (B2.isZeroAt(r) == B.isZeroAt(r));
      assert (B2.valuesAt(r) == B.valuesAt(r));
            if B2.isZeroAt(r) then 
          (             
          assert (B2.jacobianAt(r) == B.jacobianAt(r));
          ) else ();
     ));

assert  B2.isZeroAt(line)


assert(sub(B.jacobian,line)== sub(jacobian I,line))
assert (B.jacobianAt(line) == dropDegreeInfo( sub(jacobian I,line)) )


///


TEST ///
    --debug BlackBoxIdeals
    
    bbRankM = blackBoxParameterSpace( 5 ,ZZ )
    assert(bbRankM.numVariables==5);
    assert(bbRankM.coefficientRing===ZZ);

    -- assert( bbRankM.numGenerators() === null)

    rankMat := ( blackBox, point )->5 

    assert(not  bbRankM.hasPointProperty("rankMat") );

    bbRankM.registerPointProperty( "rankMat", rankMat )

    assert(  bbRankM.hasPointProperty("rankMat") );

    point  = matrix {{1,2,3,4,5}};

    assert( rankMat(bbRankM, point) == (bbRankM.pointProperty("rankMat"))(point) );

    rankMatNew := (blackBox, point)->3
    bbRankM.updatePointProperty("rankMat",rankMatNew)
    assert( rankMatNew(bbRankM,point) == (bbRankM.pointProperty("rankMat"))(point) );
    assert( rankMatNew(bbRankM,point) == (bbRankM.pointProperty(getGlobalSymbol "rankMat"))(point) );


    bbRankMNew = rebuildBlackBox bbRankM;

    assert bbRankMNew#?(global rankMat);
    assert bbRankMNew#?("rankMat");

    assert( rankMatNew(bbRankMNew,point) == bbRankMNew.rankMat(point) );

    assert(bbRankMNew.coefficientRing===ZZ);

    rankMatNew := (blackBox, point)->4 --also influences bbRankM; because rebuild does not copy; it just exports new registered properties.

    bbRankMNew.updatePointProperty("rankMat",rankMatNew)

    (bbRankM.pointProperty("rankMat"))(point);

    assert( rankMatNew(bbRankMNew,point) == bbRankMNew.rankMat(point) );
    assert( rankMatNew(bbRankMNew,point) == bbRankMNew#"rankMat"(point) );
    assert( rankMatNew(bbRankMNew,point) == (bbRankMNew.pointProperty("rankMat"))(point) );
    assert( rankMatNew(bbRankMNew,point) == (bbRankMNew.pointProperty(getGlobalSymbol "rankMat"))(point) );
    keys bbRankM
    
    valuesAt := ( blackBox, point )-> matrix {{1,2}};

    bbRankM.registerPointProperty( "valuesAt", valuesAt );
    -- that is not good; registering a point property requires a rebuild.
    bbRankM = rebuildBlackBox bbRankM
     
    assert(  bbRankM.hasPointProperty("isZeroAt") );

    assert(  bbRankM.hasPointProperty("jacobianAt") );
 

    assert( bbRankM.numGenerators() =!= null)

    assert( bbRankM.numGenerators() === 2 )
   -- bbRankMNew.numGenerators()
    
    illegalPoint := matrix {{1,2,3,4,5,6}}; 
   
    try ( bbRankMNew.rankMat(illegalPoint) ) then ( assert(false) ) else ();


    bbRankM = blackBoxParameterSpace( 5 ,ZZ/7 )

    valuesAt := ( blackBox, point )-> matrix {{1,2}};

    --bbRankM.registerPointProperty("valuesAt",valuesAt);
    bbRankM.rpp("valuesAt",valuesAt);

    bbRankM = rebuildBlackBox bbRankM;

    point  = sub(point,ZZ/7); 

    bbRankM.valuesAt(point)

    illegalPoint := sub(point,ZZ/2); 

    try ( bbRankM.valuesAt(illegalPoint) ) then ( assert(false) ) else ();

   
 ///


-- how to document 
-- blackBoxParameterSpace(ZZ,Ring) ?

doc ///
   Key
        blackBoxParameterSpace
   Headline
        create a BlackBoxParameterSpace  describing an parameter space
   Usage   
        blackBoxParameterSpace(rng)
        blackBoxParameterSpace( variableNumber, coeffRing)
   Inputs  
        rng:Ring
             a polynomial ring 
        variableNumber:ZZ
             number of parameter variables
        coeffRing:Ring
             ring to which the parameter variables belong to; 

   Outputs
        : BlackBoxParameterSpace
   Description  
        Text    
           Creates a blackbox describing an parameter space    \break            
           \break  Example:  create an {\tt BlackBoxParameterSpace } object :
        Example     
            coeffRing := ZZ/3;
            numVariables := 5;
            bbParamSpace = blackBoxParameterSpace( numVariables , coeffRing );
            rankMat := (blackBox, point)->5;
            bbParamSpace.knownMethods()
            bbParamSpace = bbParamSpace.registerPointProperty("rankJacobianAt", rankMat );
        Text
            \break Now access some propeties via the black box interface:
        Example          
            -- keys bbParamSpace
            bbParamSpace.knownAttributes()
            bbParamSpace.coefficientRing
            bbParamSpace.knownPointProperties()
            point := matrix { {1,2,3,4,5} };
            point = sub( point, coeffRing); 
            bbParamSpace.rankJacobianAt(point)
            rankMatNew := (blackBox, point)->6;
            bbParamSpace = bbParamSpace.updatePointProperty("rankJacobianAt", rankMatNew );
            bbParamSpace.rankJacobianAt(point)
///


doc ///
    Key
        BlackBoxParameterSpace
        (NewFromMethod,BlackBoxParameterSpace,Ring)
        (NewFromMethod,BlackBoxParameterSpace,Thing)
    Headline
          black boxes for parameter spaces
    Description
        Text
            A black box parameter space is used to implement parameter spaces
            with their universal families in a pointwise fashion.
            
            @TO "Singularities of cubic surfaces" @
      
            Implements an unified interface for some explicit and implicit given ideals. \break  \break
            see also @TO BlackBoxIdeal @
            \break             \break
            The simplest parameter space black box  implements  \break \break 
            Attributes:\break 
            \,\, \bullet \,{\tt numVariables }: number of variables in the parameter space. \break
            \,\, \bullet \,{\tt coefficientRing }: the ring the parameter variables belong to (or embeddable to) \break
            
            Methods:\break 
            \,\, \bullet \,{\tt knownPointProperties() }: returns a list of all known properties at a point for a blackbox \break
            \,\, \bullet \,{\tt knownMethods() }: returns a list of all known methods of the blackbox \break 
            \,\, \bullet \,{\tt knownAttributes() }: returns a list of all known attributes  \break 
            \,\, \bullet  \,{\tt registerPointProperty(propertyName, propertyMethod) }: \break 
            \,\, \, \, register a new point property for a BlackBox. \break
            \,\, \, \, e.g. evaluation at a point. \break
            \,\, \, \, expected Interface of {\tt propertyMethod} is: \break  
            \,\, \, \, \, (  {\tt blackBox }: @TO BlackBoxParameterSpace @,  {\tt point }: @TO Matrix@ )  \break 
            \break
            \,\, \, \, There are several special property names: {\tt valuesAt, jacobianAt}. \ break 
            \,\, \, \,  If one of this properties is registered or updated, it triggers update of dependent properties \break
            \,\, \, \,  e.g. registering evaluation {\tt 'valuesAt'}  will implicitly \break 
            \,\, \, \, construct  {\tt 'isZeroAt', 'numGenerators',  'jacobianAt', 'rankJacobianAt' } \break 
            \,\, \, \,  registering evaluation {\tt 'jacobianAt'}  will implicitly \break 
            \,\, \, \, construct  {\tt  'rankJacobianAt' } \break 
            \break
            \,\, \bullet  \,{\tt updatePointProperty(propertyName, propertyMethod) }: update an existing point property for a BlackBox. \break
            \,\, \bullet  \,{\tt pointProperty(propertyName) }: get a point property by name . \break 
            \,\, \bullet  \,{\tt getUpdatedBlackBox() }: return a blackbox {\tt bb } for which all recently registered properties \break 
            \,\, \, \, are accessible using the point operator: e.g. if {\tt "bettiAt"} was registered, then it is then accessible as {\tt bb.bettiAt }  \break \break  

            In addition, if the black box was created from an ideal () or an evaluation (), it has additional methods and attributes,
            see @TO BlackBoxIdeal @ 

            The blackbox interface allows implementation of some algorithms e.g. padic lifting. 
        Text
           \break  For usage example see @TO blackBoxParameterSpace@
///
            
doc ///
    Key
        "Singularities of cubic surfaces"
    Headline
        use a blackBoxParameterSpace to study the space of cubic surfaces
    Description
        Text
            A black box parameter space is used to implement parameter spaces
            with their universal families in a pointwise fashion.
            
            Let us build the parameter space of cubic surfaces with a view of 
            studying its stratification with respect to singularity type.
            
            We work in charateristic 7.            
        Example    
            K = ZZ/7
        Text
            The coordinate ring of IP^3
        Example
            R = K[x,y,z,w]
        Text  
            Make an empty blackbox which will later contain our 
            describtion of the parameter space of cubic surfaces.
            It will depend on 20 parameters, since acubic polynomial 
            in 4 variables has 20 coefficients.
        Example
            bbC = blackBoxParameterSpace(20,K);
            bbC.knownPointProperties()
        Text
            We now build the cubics from the coefficents, i.e. we
            construct the member of the universal familiy over 
            a given parameter point:
        Example
            mons3 = matrix entries transpose super basis(3,R)
            cubicAt = (point) -> matrix entries (point*mons3)
        Text 
            register this function in the black box 
        Example
            bbC = bbC.registerPointProperty("cubicAt",cubicAt);
            bbC.knownPointProperties()
        Text
            Lets test this functionality with some special cubics.
            The first example is the cubic cone. It is singular
            at (0:0:0:1):                   
        Example    
            cubicCone = matrix{{x^3+y^3+z^3}}
            coeffCubicCone = contract(transpose mons3,cubicCone)
            bbC.cubicAt(coeffCubicCone)
        Text
            The second example is the Fermat cubic. It is smooth everywhere    
        Example
            cubicFermat = matrix{{x^3+y^3+z^3+w^3}}
            coeffCubicFermat = contract(transpose mons3,cubicFermat)
            bbC.cubicAt(coeffCubicFermat)
        Text
            Now we want to implement the stratification by singularity type.
            For this we first determine the singular locus of a cubic surface:
        Example
            singularLocusAt = (bb,point) -> ideal jacobian bb.cubicAt(point)
            bbC = bbC.rpp("singularLocusAt",singularLocusAt);
            bbC.knownPointProperties()
            bbC.singularLocusAt(coeffCubicCone)   
            bbC.singularLocusAt(coeffCubicFermat)
        Text
            As a first approximation of the singularity type we use
            the degree of the singular locus
        Example
            degreeSingularLocusAt = (bb,point) -> (
                 s := bb.singularLocusAt(point);
                 if dim s == 0 then return 0;                
                 if dim s == 1 then return degree s;                 
                 if dim s >= 2 then return infinity;
              )
            bbC = bbC.rpp("degreeSingularLocusAt",degreeSingularLocusAt);
            bbC.knownPointProperties()
        Text
            Calculate the degree of the singular locus for our examples
        Example
            bbC.degreeSingularLocusAt(coeffCubicCone)
            bbC.degreeSingularLocusAt(coeffCubicFermat)
        Text
            Now the BlackBoxParameterspace hat a number of point properties
        Example
            bbC.knownPointProperties()
        Text
            These properties can now be used in a finite field experiment
            that studies the statification of our parameter space. Here is a
            simple minded version of such an experiment:
        Example
            tally apply(100,i->bbC.degreeSingularLocusAt(random(K^1,K^20))) 
        Text
            We see that there is an open stratum of smooth cubics. The
            largest closed stratum consists of those cubics with a A1 singularity.
            The package finiteFieldExperiments helps to do the bookkeeping 
            for such experiments and also provides more detailed interpretation
            of the results.

///

doc ///
    Key
        BlackBoxIdeals
    Headline
          black boxes for explicit and implicitly given ideals
    Description
        Text
            Implements an unified interface for some explicit and implicit given ideals \break  
            \break
            Purpose: use a unified interface for finite field experiments; see FiniteFieldExperiments and padicLift package
            \break
            Currently three BlackBox constructors are available:\break
            \,\,  \bullet \,   @TO blackBoxParameterSpace@  \break
            \,\,  \bullet \,   @TO blackBoxIdeal@  \break
            \,\,  \bullet \,  @TO blackBoxIdealFromEvaluation@  \break \break

            All black boxes implement the interface of @TO BlackBoxParameterSpace @ \break
            If the black box was created from an ideal or an evaluation, it has additional properties, see  @TO BlackBoxIdeal @ \break  
                        
      
    Caveat
         the black box properties are write-protected to prevent accidental modification  by the user; \break 
         however implementing write-protection leads to undesired code compexity and at the same time it is still possible to overcome the  protection. \break
         Currently adding properties to the blackBox with more than one parameter (point) is not implemented (e.g. jet computation ). \break
         Also not done yet is the implementation of the {\tt Observable } interface for the various black boxes ({\tt FiniteFieldExperiment } which will be  {\tt Observers }) \break
         Finally, the package is probably not threadsafe.
         
         
///

TEST  /// 
 debug BlackBoxIdeals
 idealBlackBoxesProtect()

 -- bblog is not defined... why ?

 BlackBoxLogger.debug("test update valuesAt property ")
 rng := ZZ/7[x]

  I = ideal(6*x)

  bb=blackBoxIdeal I

  result :=  matrix{{5}};

  point := matrix{{0_rng}};

  bb.updatePointProperty ("valuesAt", (blackBox,point)->result )  --fails

 

  assert (result==bb.valuesAt(point) );
  
 result =  matrix{{0}};

 bb.updatePointProperty ("valuesAt", (blackBox,point)->result )  --fails

  assert (result==bb.valuesAt(point) );

 
 
///

doc ///
    Key
        getEpsRing
        (getEpsRing, Ring, ZZ)
        eps
    Headline
        get a ring for jet calculations
    Usage   
        getEpsRing(R,d)
    Inputs  
        R: Ring
             the coefficient Ring
        d: ZZ 
         the precision of the Ring for jet calculations
    Outputs
        : Ring
             R[eps]/eps^{d+1}
    Description
        Text
           The advantage of using this function is that for each
           precision and coefficient ring R[eps]/eps^{d+1} is
           created only once.
        Example          
          E2 = getEpsRing(QQ,2)
          eps^2+eps^3
          E3 = getEpsRing(QQ,3)
          eps^2+eps^3
        Text
          Be aware that there are now several eps-es. This can
          lead to unexpected behavior:
        Example
          use E2; eps
        Text  
          Notice that this still gives the eps in E3.
          The underscore notation also does not help in this situation:
        Example
          try eps_E2 then "works" else "error"
        Text      
          The following works, but is not recomended since the 
          code does not check whether eps is a variable in a ring:
        Example
          sub(eps,E2)
          eps = 1
          sub(eps,E2)
        Text  
          We recommend the following:
        Example
          E2.eps
          E3.eps
          E2.eps^3
        Text
          This was implemented by hand for these rings. It does not work
          in general. 
///

doc ///
    Key
        isProbablySmoothAt
        isCertainlySingularAt
    Headline
        heuristic test of smoothness
    Usage   
        bb.isProbablySmoothAt(point)
        bb.isCertainlySingularAt(point)
    Inputs  
        bb: BlackBoxIdeal
             an black box ideal
        point: Matrix 
         the coordinates of a point
    Outputs
        : Boolean
    Description
        Text
          Checks for smoothness of a point on the
          vanishing set of the black box ideal, by
          trying to find a jet starting at the point.
          If the point is smooth on the vanishing
          set of the black box ideal, arbitray jets
          can always be found. If one or more jets are found
          the point is probably smooth. If the search
          for a jet failes only once, the vanishing
          set is certainly singular at this point.
          
          Consinder for example the cuspidal cubic:
        Example
          R = QQ[x,y]
          bbI = blackBoxIdeal ideal(x^2-y^3);
        Text
          The cuspidal cubic is singular at the origin:
        Example    
          origin = matrix{{0,0_QQ}}
          bbI.isCertainlySingularAt(origin)
        Text
          Consider a point on the cuspidal cubic different from the origin:
        Example
          otherPoint = matrix{{8,4_QQ}}
        Text
          Check whether the other point lies on the cuspidal cubic:
        Example  
          bbI.isZeroAt(otherPoint)
        Text
          We now check for smoothness:
        Example
          bbI.isCertainlySingularAt(otherPoint)
          bbI.isProbablySmoothAt(otherPoint)
        Text
          If the point is not on the vanishing set defined by
          the black box ideal, "is certainly singular" is always
          returned. (This is usefull when this property is watched
          in a finite field experiment.)
        Example
          pointNotOnCurve = matrix{{4,5_QQ}}
          bbI.isZeroAt(pointNotOnCurve)
          bbI.isCertainlySingularAt(pointNotOnCurve)
          jetAt(bbI,pointNotOnCurve,1,1)
      
///

doc ///
    Key
        jetAt
        (jetAt, BlackBoxParameterSpace, Matrix, ZZ, ZZ)
    Headline
        find jets on the varieties defined by a black box
    Usage   
        jetAt(bb,point,prec,n)
    Inputs  
        bb: BlackBoxIdeal
             an black box ideal
        point: Matrix 
             the coordinates of a point
        prec: ZZ
             the precision of the desired jet
        n: ZZ
             the number of trial used to find a jet    
    Outputs
        : MutableHashTable
    Description
        Text
          Tries to find a jet starting at a given point on 
          a variety given by a black box.
          
          Consinder for example a nodal cubic:
        Example
          Fp = ZZ/101
          R = Fp[x,y]
          I = ideal(x^2-y^2+x^3)
          bbI = blackBoxIdeal I;
        Text
          Consider a point on the nodal cubic different from the origin:
        Example
          point = matrix{{3,6_Fp}}
        Text
          Check whether the other point lies on the nodal cubic:
        Example  
          bbI.isZeroAt(point)
        Text
          We now look for a jet:
        Example
          j = jetAt(bbI,point,3,1)
        Text
          The defining equations of the ideal indeed vanish on the jet:
        Example
          sub(I,j#"jet")
        Text
          At the origin the nodal cubic is singular. Short jets can be found,
          but not long ones:
        Example
          origin = matrix{{0,0_Fp}}
          jetAt(bbI,origin,1,1)  
          jetAt(bbI,origin,2,1)  
          jetAt(bbI,origin,3,1)  
        Text
          Notice that the search for fails at length 2 most of the time,
          since the singularity has multiplicity 2. If one tries
          long enough a longer jet can be found (lying on one
          of the branches at the origin):
        Example        
            jetAt(bbI,origin,3,200)  
      
///


doc ///
    Key
        rebuildBlackBox
    Headline
        provide new added black box properties via point('.') operator
    Usage   
        rebuildBlackBox(bb)
    Inputs  
        bb: BlackBoxIdeal
             a black box ideal
        bb: BlackBoxParameterSpace
             a black box parameter space
    Outputs
        : BlackBoxIdeal
        : BlackBoxParameterSpace
    Description
        Text
          Provide via registerPointProperty() added new black box properties
          via point('.') operator (syntactic sugar)  
        Example
            K = ZZ/101
            bbC = blackBoxParameterSpace(2,K);

            valueAt = (point) -> ( 
                 R = K[x,y];
                 I = ideal(x^2-y^2+x^3);
                 return sub(I,sub(point,R) )
             )
        Text
            register this function in the BlackBox
        Example
            updatedbbC = bbC.registerPointProperty("valueAt",valueAt);
        Text
            in case we do not use updatedbbC we cannot access 'valueAt' via point operator:
        Example
            point = matrix{{3,6_K}};
            try bbC.valueAt(point) else "error"
            updatedbbC.valueAt(point) -- ok
        Text
            however, we can get construct another updated black box using 'rebuildBlackBox()'
        Example
            bbC = rebuildBlackBox(bbC);
            bbC.valueAt(point) -- ok
///

doc ///
    Key
        "isZeroAt"
    Headline
         Check if a given point lies on the vanishing set defined by a black box ideal.
    Usage   
         bbI.isZeroAt(point)
    Inputs  
        bbI: BlackBoxIdeal
        point: Matrix
             coordinates of a point
    Outputs
        : Boolean
    Description
        Text
          This is a point property.
          Checks if the the point lies on the vanishing set defined by the
          black box ideal. 
        Example
          R = QQ[x,y]
          bbI = blackBoxIdeal ideal(x^2-y^3);
          bbI.isZeroAt(matrix{{0,0_QQ}})
          bbI.isZeroAt(matrix{{8,4_QQ}})            
          bbI.isZeroAt(matrix{{8,5_QQ}})
///

doc ///
    Key
        "valuesAt"
    Headline
        evaluate the generators of black box ideal at a given point 
    Usage   
        bbI.valuesAt(point)
    Inputs  
        bbI: BlackBoxIdeal
        point: Matrix
             coordinates of a point
    Outputs
        : Boolean
    Description
        Text
          This is a point property. It evaluates the generators of the
          black box ideal at the point
        Example
          R = QQ[x,y]
          bbI = blackBoxIdeal ideal(x^2,y^3);
          bbI.valuesAt(matrix{{2,3_QQ}})
        Text
          valuesAt is a special pointProperty
          in the sense that the package uses them to compute other
          point properties. Therefore
          updating valuesAt triggers updates for 
          isZeroAt, numGenerators, jacobianAt, rankJacobianAt. 
    Caveat
          This works only with black box ideals, since they contain an algorithm
          that can evaluate the generators of the black box ideal. A black box parameter spaces
          might contain only an algorithm that checks whether all generators vanish. 
          This happens for example if one considers the moduli space of singular cubics. 
          One can check whether a
          given cubic is singular without calculating the value of the 
          corresponding discriminant. 
///

doc ///
    Key
        "jacobianAt"
    Headline
        evaluate the jacobian matrix of a black box ideal at a given point 
    Usage   
        bbI.jacobianAt(point)
    Inputs  
        bbI: BlackBoxIdeal
        point: Matrix
             coordinates of a point
    Outputs
        : Boolean
    Description
        Text
          This is a point property. It evaluates the jacobian matrix of the
          black box ideal at the point
        Example
          R = QQ[x,y]
          I = ideal(x^2-y^3);
          bbI = blackBoxIdeal I;
          point = matrix{{8,4_QQ}}
          bbI.isZeroAt(point)
          bbI.jacobianAt(point)
          sub(jacobian I,point)
          bbI.isProbablySmoothAt(point)
        Text
          The cuspidal cubic considered above is singular at the
          origin. Therefore the jacobian matrix vanishes there:
        Example
          origin = matrix{{0,0_QQ}}
          bbI.isZeroAt(origin)
          bbI.jacobianAt(origin)
          bbI.isCertainlySingularAt(origin)
        Text
          jacobianAt is a special pointProperty
          in the sense that the package uses it to compute other
          pointProperties. Therefore an update of 
          jacobianAt triggers an update 
          for 'rankJacobianAt'.  
    Caveat
          This works only with black box ideals, since they contain an algorithm
          that can evaluate the generators of the black box ideal. A black box parameter spaces
          might contain only an algorithm that checks whether all generators vanish. 
          This happens for example if one considers the moduli space of singular cubics. 
          One can check whether a
          given cubic is singular without calculating the value of the 
          corresponding discriminant. 
///

doc ///
    Key
        "registerPointProperty"
        "rpp"
    Headline
        register a new point property in a black box.
    Usage   
        bbI = bbI.registerPointProperty(name,propertyAt)
        bbI = bbI.rpp(name,propertyAt)
    Inputs  
        bbI: BlackBoxParameterSpace
        name : String
          name of the new point property
        propertAt: Function 
          that takes coordinates of a point and returns anything.
    Outputs
        : BlackBoxParameterSpace
    Description
        Text
          rpp and 
          registerPointProperty are synonymous. rpp is provided
          to save typing time...
          
          This method is used to register new property in a
          blackBoxIdeal or a blackBoxParameterSpace. 
          
          Lets for example build a black box parameter space
          for cubic surfaces in IP^3 over the finite field with
          7 Elements.
        Example
          Fp = ZZ/7
          R = Fp[x,y,w,z]
        Text
          there are 20 monomials of degree 3 in 4 variables
        Example
          mons3 = basis(3,R)
        Text
          Therefore we need a 20 dimensional parameter space
        Example
          bbC = blackBoxParameterSpace(20,Fp);
        Text
          This has no known point properties 
        Example
           bbC.knownPointProperties()
        Text
          We now make a function that constructs a cubic
          polynomial from 20 parameters
        Example
          cubicAt = (point) -> point*transpose mons3
          cubicAt(matrix{{1_Fp,19:0}})
          fermatPoint = matrix {{1_Fp,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1}}
          cubicAt(fermatPoint)
        Text
          To have this available inside the blackbox we need to register
          it
        Example
          bbC = bbC.registerPointProperty("cubicAt",cubicAt);
          bbC.knownPointProperties()
        Text
          Now we can use the property from the black box
        Example
          bbC#"cubicAt"(fermatPoint)
          bbC.cubicAt(fermatPoint)
        Text
          Registering a point property is useful when the black box
          is used in a finite field experiment. Registered
          point properties are available to an experiment while it is
          running. Also it is useful to register for bookkeeping reasons,
          so that for example knownPointProperties will return
          the correct answer.
///

doc ///
    Key
        "updatePointProperty"
        "upp"
    Headline
        change a point property in a black box.
    Usage   
        bbI.updatePointProperty(name,propertyAt)
        bbI.upp(name,propertyAt)
    Inputs  
        bbI: BlackBoxParameterSpace
        name : String
          name of point property to be changed
        propertAt: Function 
          that takes coordinates of a point and returns anything.
    Outputs
        : BlackBoxParameterSpace
    Description
        Text
          upp and 
          updatePointProperty are synonymous.
          
          This method is used to change a property of 
          blackBoxIdeal or a blackBoxParameterSpace. The need
          for this arises usually while programming when one
          realizes that a programming mistake was made.  
          
          Lets look at a stupid, but illustrative example:       
        Example
          bbC = blackBoxParameterSpace(2,QQ);
          bbC = bbC.rpp("product",(point) -> sum flatten entries point);
          bbC.product(matrix{{5,6_QQ}})
          bbC.upp("product",(point) -> product flatten entries point);
          bbC.product(matrix{{5,6_QQ}})
        Text
          It is also possible to update pointProperties that are
          predefined by the package. Here
          valuesAt and jacobianAt are special pointProperties,
          in the sense that the package uses them to compute other
          pointProperties. Therefore
          updating valuesAt triggers updates for 
          isZeroAt, numGenerators, jacobianAt, rankJacobianAt.
          Similarily updating jacobianAt triggers an update 
          for 'rankJacobianAt'.  
    SeeAlso
         registerPointProperty
         rpp
         valuesAt
         jacobianAt
///

doc ///
    Key
        "hasPointProperty"
    Headline
        check whether a point property of a black box is defined
    Usage   
        bbI.hasPointProperty(name)
    Inputs  
        bbI: BlackBoxParameterSpace
        name : String
          name of point property to be checked
    Outputs
        : Boolean
    Description
        Text
          check whether a point property is defined.
          
          Every BlackBoxIdeal has the property "jacobianAt" since
          it has (at least implicitly) access to a representation of 
          the generators of the ideal:
        Example
          R = QQ[x,y]
          bbI = blackBoxIdeal ideal (x^2-y^3);
          bbI.hasPointProperty("jacobianAt")
          bbI.knownPointProperties()
        Text
          A blackBoxParameterSpace usually does not have 
          the property "jacobianAt" since such a parameter space
          does not even have an implicit representation of 
          the equations.
        Example
          bbParam = blackBoxParameterSpace(2,QQ);     
          bbParam.hasPointProperty("jacobianAt") 
          bbParam.knownPointProperties()
        Text
          To illustrate why this can happen think for example of 
          the space of
          singular cubics. In a BlackBoxParameterSpace one
          would simply test the smoothness of a give cubic
          via Groebner basis calculation. This does not automatically
          give rise to a representation of the corresponding
          Diskriminant. 
    SeeAlso
         knownPointProperties
         registerPointProperty
///

doc ///
    Key
        "knownAttributes"
        "numVariables"
    Headline
        list the attributes of a black box 
    Usage   
        bbI.knownAttributes()
    Inputs  
        bbI: BlackBoxParameterSpace
    Outputs
        : List       
    Description
        Text
          Every BlackBoxIdeal has some attributes
          provided by the package 
          (no new attributes can be defined by the user).
          The difference between attributes and methods is,
          that an attribute is a constant while a property is 
          a function.
        Example
          R = QQ[x,y]
          bbI = blackBoxIdeal ideal (x^2-y^3);
          bbI.knownAttributes()
        Text
          Lets look at these attributes in turn
        Example
          bbI.coefficientRing
          bbI.equations
          bbI.ideal
          bbI.jacobian
          bbI.numVariables
          bbI.ring
          bbI.type
        Text
          The type can also be "BlackBoxParameterSpace"  
        Example
          bbI.unknowns
        Text
          Lets now look at a blackBoxIdeal defined by an
          evaluation. The standart example is the determinant
          of a matrix:
        Example
          M = (point) -> (point_{0,1,2}||point_{3,4,5}||point_{6,7,8})
          phonePoint = matrix{{1,2,3,4,5,6,7,8,9_QQ}}
          M(phonePoint)
          detM = (point) -> matrix{{det(M(point))}}
          detM(phonePoint)   
          S = QQ[m_1..m_9]  
          bbE = blackBoxIdealFromEvaluation( S, detM );
          bbE.valuesAt(phonePoint)
          bbE.knownAttributes()
        Text
          Notice that "equations" and "jacobian" are missing, since
          no explicit equations of the blackBoxIdeal are provided.
        Example
          bbP = blackBoxParameterSpace(2,QQ);
          bbP.knownAttributes()
        Text
          For a blackPointParameterSpace "ring" and "unknowns" are 
          missing since there are no equations (not even implicit ones)  
    SeeAlso
         knownPointProperties
         knownMethods
///

doc ///
    Key
        "knownMethods"
    Headline
        list the methods of a black box 
    Usage   
        bbI.knownMethods()
    Inputs  
        bbI: BlackBoxParameterSpace
    Outputs
        : List       
    Description
        Text
          Every BlackBoxIdeal has some methods
          provided by the package 
          (no new attributes can be defined by the user).
          The difference between attributes and methods is,
          that an attribute is a constant while a property is 
          a function.
        Example
          R = QQ[x,y]
          bbI = blackBoxIdeal ideal (x^2-y^3);
          bbI.knownMethods()
        Text
          Each of these methods has their own documentation node.
    SeeAlso
         hasPointProperty
         knownAttributes
         knownPointProperties
         knownPointPropertiesAsSymbols
         registerPointProperty
         rpp
         setSingularityTestOptions
         updatePointProperty
         upp
///

doc ///
    Key
        "knownPointProperties"
    Headline
        list the pointProperties of a black box 
    Usage   
        bbI.knownPointProperties()
    Inputs  
        bbI: BlackBoxParameterSpace
    Outputs
        : List       
    Description
        Text
          A pointProperty of a black box is a function that depends
          only on the coordinates of a point in the parameter space. 
        Example
          R = QQ[x,y]
          bbI = blackBoxIdeal ideal (x^2-y^3);
          bbI.knownPointProperties()
          bbI.isZeroAt(matrix{{0,0_QQ}})
          bbI.jacobianAt(matrix{{1,1_QQ}})
          bbI.isProbablySmoothAt(matrix{{1,1_QQ}})
          bbI.isCertainlySingularAt(matrix{{0,0_QQ}})
        Text
          Each of these pointProperties has their own documentation node.
          New pointProperties can be defined by user using registerPointProperty.
          Exisiting pointProperties can be changed by using 
          updatePointProperty. 
          
          valuesAt and jacobianAt are special pointProperties,
          in the sense that the package uses them to compute other
          point properties. Therefore
          updating valuesAt triggers updates for 
          isZeroAt, numGenerators, jacobianAt, rankJacobianAt.
          Similarily updating jacobianAt triggers an update 
          for 'rankJacobianAt'.
    SeeAlso
      isCertainlySingularAt
      isProbablySmoothAt
      isZeroAt
      jacobianAt
      rankJacobianAt
      valuesAt
      registerPointProperty
      updatePointProperty
 ///



doc ///
    Key
        "rankJacobianAt"
    Headline
        determine the rank of the jacobian matrix of a black box ideal at a given point 
    Usage   
        bbI.rankJacobianAt(point)
    Inputs  
        bbI: BlackBoxIdeal
        point: Matrix
             coordinates of a point
    Outputs
        : Boolean
    Description
        Text
          This is a point property. It evaluates the jacobian matrix of the
          black box ideal at the point and determines its rank.
        Example
          R = QQ[x,y]
          I = ideal(x^2-y^3);
          bbI = blackBoxIdeal I;
          point = matrix{{8,4_QQ}}
          bbI.isZeroAt(point)
          bbI.jacobianAt(point)
          bbI.rankJacobianAt(point)
          bbI.isProbablySmoothAt(point)
        Text
          The cuspidal cubic considered above is singular at the
          origin. Therefore the jacobian matrix hat rank 0 there:
        Example
          origin = matrix{{0,0_QQ}}
          bbI.isZeroAt(origin)
          bbI.jacobianAt(origin)
          bbI.rankJacobianAt(origin)
          bbI.isCertainlySingularAt(origin)
        Text
          This point property is usefull when running experiments on
          black boxes that have serveral components of unknown 
          codimension. The rank of the jacobian matrix gives a upper bound on 
          the codimension of the components containing the point (with
          equality if the vanishing set is smooth at the point). Sorting 
          the number of points found in an experiment by the rank of their
          jacobian matrices
          helps to estimate the number of components of the vanishing set
          of the black box in each codimension.   
    Caveat
          This works only with black box ideals, since they contain an algorithm
          that can evaluate the generators of the black box ideal. A black box parameter spaces
          might contain only an algorithm that checks whether all generators vanish. 
          This happens for example if one considers the moduli space of singular cubics. 
          One can check whether a
          given cubic is singular without calculating the value of the 
          corresponding discriminant. 
///

end
-- need test for randomIterator (for fixed error : to less trials  )
--@TO2{jetAt,"jet"}@ computations at a point independently of the ideal representation. \break \break 


--        (clearCoeffDenominators, Ideal )    


