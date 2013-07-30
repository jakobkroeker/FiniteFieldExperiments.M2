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
     Headline => "black boxes for explicit and implicitly given ideals",
     DebuggingMode => true
)

export {
    --clearCoeffDenominators,
    "clearCoeffDenominators",
    blackBoxIdeal,
    blackBoxIdealFromEvaluation,
    blackBoxIdealFromProperties,
    getEpsRing,
    JetAt,
    createLogger
}

--undocumented {
--
--}

jacobianAt := global jacobianAt;
jacobianAt = global jacobianAt;

idealBlackBoxesProtect = ()->
(

protect jacobianAt;
protect jetAt;
protect jacobianMatrix;
protect transposedJacobianAt;
protect transposedJacobian;
protect valuesAt;
protect unknownIsValid;
protect createBlackBoxIdeal;
protect getEquations;
protect numVariables;
protect imageRank;
protect isZeroAt;
protect pointProperties;
protect setRing;
protect registerPointProperty;
protect registerAnonymousPointProperty;
protect setPointProperty;
protect setValuesAt;
protect setImageRank;
protect checkInputPoint;
protect deduceImageRank;
protect propertiesAt;
protect setPropertiesAt;
protect setIsZeroAt;
protect dropDegreeInfo;
protect bareJacobianAt;
protect rebuild;
protect setThis;
protect clearInternal;
protect unknowns;
protect pointProperty;
protect updatePointProperty;
protect setJacobianAt;
protect equations;
protect setLevel;
protect pointPropertyByName;
protect hasPointProperty;

protect knownPointProperties;
);

--todo: fix dublicate code,  -  padicLiftProtect and padicLiftExport
idealBlackBoxesExport = ()->
(
    exportMutable( eps );
    exportMutable( jacobianAt);
    exportMutable( jetAt);

    exportMutable( jacobianMatrix);
    exportMutable( transposedJacobian);
    exportMutable( transposedJacobianAt);
    exportMutable( valuesAt);
    exportMutable( unknownIsValid);    
    exportMutable( createBlackBoxIdeal);
    exportMutable( getEquations);                
    exportMutable( numVariables);  
    exportMutable( imageRank);
    exportMutable( isZeroAt);      
    exportMutable( pointProperties);  
    exportMutable( setRing);  
    exportMutable( registerPointProperty); 
    exportMutable(  registerAnonymousPointProperty );
    exportMutable( setPointProperty);
    exportMutable( setValuesAt);    
    exportMutable( setImageRank );
    exportMutable(  checkInputPoint);
    exportMutable( deduceImageRank );
    exportMutable( propertiesAt );
    exportMutable( setPropertiesAt );
    exportMutable( setIsZeroAt );
    exportMutable( dropDegreeInfo );
    exportMutable( bareJacobianAt );
   exportMutable( rebuild );
   exportMutable( setThis );
   exportMutable( clearInternal );
 exportMutable( unknowns );
 exportMutable( pointProperty );
 exportMutable( updatePointProperty );
 exportMutable( setJacobianAt );
 exportMutable( equations );
 exportMutable( setLevel );
 exportMutable( pointPropertyByName );
 exportMutable( hasPointProperty );
  exportMutable( knownPointProperties );
)



String ? Symbol := (str,symb)->
(
  1 ? 2
);

Symbol ? String := (str,symb)->
(
  2 ? 1
);

-- swith between protect and export - both are not possible!

--idealBlackBoxesProtect() -- protect the symbols for checking package correctness: no symbol variables should be overwritten by accident!
idealBlackBoxesExport(); -- export the symbols to make the package work 


needsPackage "SimpleDoc";
needsPackage "Text";

-- Estring = new Type of String - goal was to distinguies string keys from symbol keys without changing string visuaulization (net())
-- but that turned out to be impossible, because e.g. operator '#' cannot be overloaded.


--net (String) := Net =>(str)->(   strl := "\""| str| "\"" ;   return stack separate strl; );

savedEpsRings := new MutableHashTable;

getEpsRing = method();

getEpsRing(Ring, ZZ) := Ring => (coeffring, epsDim)->
(
    --eps := null;
    --eps = symbol eps; -- w
    -- how to use a global symbol
    --geps := getGlobalSymbol(BlackBoxIdeals.Dictionary, "eps"); 
    geps := getSymbol "eps"; 
    if (epsDim<0) then error("expected epsDim>0 ");
    
    if not (savedEpsRings#?(coeffring,epsDim) ) then 
     (
        polRing:=coeffring[geps];
        geps=(gens(polRing))#0;
        savedEpsRings#(coeffring,epsDim) = polRing/geps^(epsDim+1);    
    ); 
    return savedEpsRings#(coeffring, epsDim);
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



createLogger = ()->
(
    logger := new MutableHashTable;
    --
    level := 0;
    logger.setLevel = method();
    --
    logger.setLevel(ZZ) := (pLevel)->
    ( 
       assert(level>=0 and level<5);
       level = pLevel;
    );
    --    
    logger.log = method();
    logger.log (ZZ,String) := (msgLvl,msg )->
    (
       if (msgLvl<=level) then
       print( "--" | msg );
    );
    --
    return new HashTable from logger;
);


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
             ideal with integer coefficients with the same zero set as the input ideal
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



deduceImageRank := (blackBox)->
(
     valuesAt:=null; 
     valuesAt = global valuesAt;

     assert( blackBox#?(global valuesAt) );
     assert( blackBox#?valuesAt );

     --assert( blackBox#?(BlackBoxIdeals.Dictionary#(toString symbol valuesAt)) );  -- ok
     
       computed  := false; 
        maxTrials := 100;
        currTrial := 0;
        rng := blackBox.coefficientRing;
        imageRank := null;
        while imageRank===null and currTrial<maxTrials do
        (
          try (
              tmppoint := matrix random(rng^1,rng^(blackBox.numVariables) );
              valuesMatrix := blackBox.valuesAt( tmppoint );
              --print valuesMatrix;
              assert (numRows valuesMatrix==1);
              imageRank = numColumns valuesMatrix;

          )         then ( computed=true ) else();
               
           
          currTrial = currTrial+1;
        );
        if not computed then error "error: failed to deduce image dim";
        return imageRank;
     );



deducedJacobianAt = (blackBox,point)->
(
  
    valuesAt := null; imageRank := null;
    valuesAt = global valuesAt;
    imageRank = global imageRank;
    assert( blackBox#?(global valuesAt) );
    assert( blackBox#?(global imageRank) ); 

    rngPoint := ring point;
    numVariables := blackBox.numVariables ;      
    if (not ( blackBox.valuesAt( point )==0))  then  error("point does not belong to the ideal ! ");
    eps := null;
    eps = symbol eps;
    epsRng := rngPoint[eps]/eps^2;
    eps = (gens epsRng)#0;

    jacobianMatrixAt := mutableMatrix( rngPoint ,  numVariables, blackBox.imageRank() );
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

dropDegreeInfo := method();
dropDegreeInfo (Matrix) := Matrix=> (mat)->
(
   return map( (ring mat)^(numRows mat) ,(ring mat)^(numColumns mat),mat);
);



JetAtSingleTrial = method();

-- todo: optionally pass jacobianAt?
-- problem No 1: need same eps for all BlackBoxes?
-- problem No 2: need same epsRng for all ? 

JetAtSingleTrial(HashTable, Matrix, ZZ) := MutableHashTable => ( blackBox,  point, liftDepth )  ->
(

    print ("liftDepth",liftDepth);

    if not (blackBox.isZeroAt(point)) then error(" function is not zero at the point ");


    liftingFailed := false;
    failedLiftDepth := null;

    jacobianM2Transposed := transpose blackBox.jacobianAt(point) ;
    
    jtColumns := numColumns jacobianM2Transposed ;
    jtRows    :=  numRows jacobianM2Transposed ;

    if not ( rank jacobianM2Transposed == min( jtRows,jtColumns ) ) then 
       error ("computing jets: point is not smooth !");

    jacobianKernel := generators kernel jacobianM2Transposed ;

    epsRng := getEpsRing( blackBox.coefficientRing,  1 );
    eps = (gens epsRng)#0;

    coeffRng := (blackBox.coefficientRing); -- klammern sind Notwendig!

     
    rnd := random( coeffRng^(numColumns(jacobianKernel)), coeffRng^1 );
    if (numColumns(jacobianKernel)>0) then 
    (   
        while  zero(rnd) do
        (
            rnd = random( coeffRng^numColumns(jacobianKernel), coeffRng^1 );
        );
    );
    
    lenghtOneLift := sub(point,epsRng) + sub( jacobianKernel*rnd, epsRng) *eps;

    lift := lenghtOneLift;

    epsPrecision := 1;     -- first lift should succeed for a smooth point! 
    lastLiftNr    := epsPrecision ; 

    for  epsPrecision in 2..liftDepth do 
    (
        values := blackBox.valuesAt( lift );

        rightHandSide := matrix mutableMatrix( coeffRng, numColumns values ,1 );
        
        -- notwendig:
   
        if not zero(values) then 
            rightHandSide = transpose last coefficients (values, Monomials=>{ eps^epsPrecision });

 

        --transposedM2JacobianAtSolution =  transpose blackBox.jacobianAt( lift ) ;

        -- try is not good for debugging
        --try (
                 x := rightHandSide // jacobianM2Transposed ;
                 --jacobianKernel = generators kernel jacobianM2Transposed ;
                 x = x + jacobianKernel* random(coeffRng^(numColumns(jacobianKernel)), coeffRng^1 );

        -- ) else (
         --    failedLiftDepth = epsPrecision;
         --    liftingFailed=true; break; 
        --);
        x = transpose x;

        epsRng = getEpsRing( coeffRng, epsPrecision  );
        eps = (gens epsRng)#0;
        
        lift = sub( x, epsRng ) *eps^epsPrecision + sub (lift, epsRng );
    );
    -- todo: create a datatype or a hashTable

    retVal := new HashTable from { "succeeded" => not liftingFailed, "lift" => lift, "failedLiftDepth" =>failedLiftDepth };
    return retVal;
)
    
JetAt = method();

JetAt(HashTable, Matrix, ZZ, ZZ) := MutableHashTable => ( blackBox,  point, liftDepth, trials )  ->
(
    
    for i in 1..trials do
    (
        jetTrialResult := JetAtSingleTrial ( blackBox,  point, liftDepth);

        if (jetTrialResult#"succeeded") then 
            return jetTrialResult;
    );
    return  new HashTable from { "succeeded" =>false, "lift" => null };
    
)


JetAtWrapper = method( Options=>{"liftDepth" =>4, "trials"=>5} );

JetAtWrapper(HashTable, Matrix)  := MutableHashTable => opts -> ( blackBox,  point )  ->
(
    return JetAt( blackBox,  point, opts#"liftDepth", opts#"trials");
);


basicBlackBox = method();
--basicBlackBox (Ring) :=  HashTable =>( parrng ) ->
basicBlackBox(ZZ,Ring) := HashTable => ( numVariables, coeffRing ) ->
(
    
   blackBox := new MutableHashTable;   
   pointProperties := new HashTable  ;
   blackBox.coefficientRing = coeffRing;
   blackBox.numVariables = numVariables;
   imageRank := null;


   blackBox.imageRank = () ->
   (  return imageRank; );
   
     

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


   

     blackBox.pointProperty=(prop)->
     (
        if not pointProperties#?prop then return null;

        return pointProperties#prop;
     );

  

   blackBox.hasPointProperty=(propertyName)->
   (
       return pointProperties#?propertyName;
   );


   setPointProperty := method();
   setPointProperty ( Symbol, Function, Function) := Thing => ( propertySymbol, propertyMethod, preconditionsTest )->
   (
      -- propertySymbol := getGlobalSymbol propertyName;
      propertyName := toString propertySymbol;

      assert(propertyName=!=null); 
      
      packageSymbol := null;
      if BlackBoxIdeals.Dictionary#?propertyName then 
        packageSymbol  = BlackBoxIdeals.Dictionary#propertyName;
    

      preconditionsTest(blackBox); -- e.g. check for method presence or test an example.

      pointProperties = new MutableHashTable from pointProperties;
 
      
          pointProperties#propertyName = (point )-> ( 
          checkInputPoint(point);
          try { return propertyMethod(  point); } then {  return propertyMethod(  point); } else {     return propertyMethod(blackBox,   point);  };
     );

     pointProperties#propertySymbol = pointProperties#propertyName;
  
     if packageSymbol=!=null then 
       pointProperties#packageSymbol = pointProperties#propertyName;
     
     if mutable blackBox then 
     (
         blackBox#propertySymbol        = (point)->( (blackBox.pointProperty(propertySymbol))(point) );
         blackBox#propertyName          = (point)->( (blackBox.pointProperty(propertyName))(point) );
         if packageSymbol=!=null then 
           blackBox#packageSymbol        = (point)->( (blackBox.pointProperty(packageSymbol))(point) );
     );

     pointProperties = new HashTable from pointProperties;
  
   );

   blackBox.knownPointProperties = ()->
   (
      return keys pointProperties;
   );


   setPointProperty( Symbol, Function) := Thing => ( propertySymbol, propertyMethod )->
   (
      setPointProperty( propertySymbol, propertyMethod, (varietyIdeal)->() );
      return null;
   );

   valuesAtWrapper := ( pValuesAt, point) ->
   (
      result := pValuesAt(point);
      assert(numRows result==1);
      return result;
   );

   setIsZeroAt := (pIsZeroAt) ->
   ( 
      --blackBox.setPointProperty( global isZeroAt , pIsZeroAt );
       setPointProperty( global isZeroAt , pIsZeroAt );
   );

   setJacobianAt := (jacobianAt) ->
   (
      blackBox.jacobianAt =jacobianAt;  
      --blackBox.setPointProperty( global jacobianAt , blackBox.jacobianAt );
      setPointProperty( global jacobianAt , blackBox.jacobianAt );
      blackBox.bareJacobianAt = (point)->
      (
         return dropDegreeInfo( blackBox.jacobianAt(point) );
      );
      setPointProperty( global bareJacobianAt , blackBox.bareJacobianAt );
   );


   setValuesAt := (pValuesAt) ->
   (      
        print "setValuesAt: updates (isZeroAt, jacobianAt)" ;      

       localValuesAt := (point)->return valuesAtWrapper(pValuesAt, point ) ;
       -- when using valuesAt instead of localValuesAt we get the wrong  (symbol valuesAt) (local valuesAt)
       --blackBox.setPointProperty( global valuesAt  ,  localValuesAt );
       setPointProperty( global valuesAt  ,  localValuesAt );
       
    
       imageRank =   deduceImageRank(blackBox)  ; --depends on valuesAt.

      print ( "updated blackBox.imageRank to " | toString blackBox.imageRank );   

       setIsZeroAt(
          (point)->( return blackBox.valuesAt(point)==0 ;) 
        );

       ----- jacobian at:
       setJacobianAt ( (point)->deducedJacobianAt( blackBox, point ) );  
   );
 

   outerSetPointProperty := ( propertySymbol, propertyMethod)->
   (
        propertyName := toString  propertySymbol;
         if propertyName==="isZeroAt" then 
              return setIsZeroAt(propertyMethod); --probably not necessary

         if propertyName==="valuesAt" then 
             return setValuesAt(propertyMethod);  -- triggers initialization of 'isZeroAt' and 'jacobianAt'
     
         if propertyName==="jacobianAt" then 
             return setJacobianAt(propertyMethod); -- triggers initialization of 'bareJacobianAt'

           setPointProperty( propertySymbol, propertyMethod );
   );

   --updating ValuesAt will not trigger updating isZeroAt and jacobianAt 
   -- 
   blackBox.updatePointProperty=( propertyName, propertyMethod )->
   (
      propertySymbol := getGlobalSymbol propertyName;

      if  (    pointProperties#?propertySymbol  and   pointProperties#?propertyName ) then
      (  
           outerSetPointProperty( propertySymbol, propertyMethod );
      ) 
      else ( error "property"| toString propertyName | "does not exist.");
   );



    -- external function 

   blackBox.registerPointProperty = method();

   -- was passiert bei Zirkelschlüssen?
   blackBox.registerPointProperty(String, Symbol, Function, Function) := Thing => ( propertyName, propertySymbol, propertyMethod, preconditionsTest )->
   (
      preconditionsTest(blackBox); -- e.g. check for method presence or test an example.
      
      if  ( not  pointProperties#?propertySymbol  and  not  pointProperties#?propertyName 
      and not  blackBox#?propertySymbol          and   not  blackBox#?propertyName  ) then 
      (
           outerSetPointProperty( propertySymbol, propertyMethod );
      )
      else error(" method "| propertyName |" is already registered");

   );


   blackBox.registerPointProperty(String, Function, Function) := Thing => ( propertyName, propertyMethod, preconditionsTest )->
   (
    
     -- propertySymbol := getGlobalSymbol propertyName; --leads to error ???
      propertySymbol := getSymbol propertyName;
     blackBox.registerPointProperty(  propertyName,propertySymbol, propertyMethod, preconditionsTest )
   );

   -- propertyMethod is expected to get two parameters, (blackBox,point)

   blackBox.registerPointProperty(String, Function) := Thing =>( propertyName, propertyMethod )->
   (
   

      emptyFkt := (blackBox)->(return null);
      blackBox.registerPointProperty( propertyName, propertyMethod, emptyFkt );
   );

   
   
   blackBox.clearInternal = ()-> 
   (
      remove( blackBox, getGlobalSymbol "setThis" );
      remove( blackBox, getSymbol "setThis" );
      remove( blackBox, symbol setThis );
      assert( not blackBox#?(symbol setThis) );
      assert( not blackBox#?( getSymbol "setThis" ) );

      remove( blackBox, getGlobalSymbol "clearInternal" );
      remove( blackBox, getSymbol "clearInternal" );
      remove( blackBox, symbol clearInternal );
      assert( not blackBox#?(symbol clearInternal) );
      assert( not blackBox#?(getSymbol "clearInternal") );
   );

   blackBox.setThis = (bb)->
   (
       blackBox = bb;
   );

   blackBox.rebuild = ()->
   (
       bb := new MutableHashTable from  blackBox;
       blackBox = bb;

       for  property in blackBox.knownPointProperties() do
       (
          if not bb#?property then 
          (
              --bb#property = pointProperties#property;
              bb#property = (point)->( blackBox.pointProperty(property)(point) );
          )
          else 
      (
        
       );         
       );

       return new HashTable from bb;
   );

   return blackBox;
)


basicBlackBox(Ring) := HashTable => ( pRing ) ->
(

   blackBox := basicBlackBox(#(gens pRing), coefficientRing pRing);

   blackBox.ring = pRing;
   blackBox.unknowns = gens blackBox.ring;
   assert( blackBox.numVariables == #blackBox.unknowns );

   blackBox.unknownIsValid = (unknown)->
   (
        if not ( blackBox.ring === ring unknown) then 
        ( 
            print "the unknown is not element of the equations ideal ring";
            return false;
        );
        return true;
   );

   --blackBox.setThis(blackBox);
   return blackBox;
)


blackBoxIdeal = method();

blackBoxIdeal(Ring) := HashTable => ( pRing ) ->
( 
    blackBox := basicBlackBox( pRing);
    blackBox.setThis(blackBox);
    blackBox.clearInternal();
    return new HashTable from blackBox;
);


--basicBlackBox(ZZ,Ring) := HashTable => ( numVariables, coeffRing ) ->
--(
--    assert ( numVariables>0 );
--    a := null;
--    a = symbol a;
--    rng := coeffRing[a_1..a_numVariables];
--    return basicBlackBox( rng );
--)



blackBoxIdeal(ZZ, Ring) := HashTable => ( numVariables, coeffRing )  ->
(
    blackBox := basicBlackBox( numVariables, coeffRing );
    blackBox.setThis(blackBox);
    blackBox.clearInternal();
    return new HashTable from blackBox;
)




-- todo: how to check, if 'ring equationsIdeal' is not a quotient ring?
blackBoxIdeal (Ideal) := HashTable =>(equationsIdeal)->
(   
     blackBox :=  basicBlackBox( ring equationsIdeal );

     blackBox.setThis(blackBox);
   
     -- maybe blackBox.addProperty( ideal, equationsIdeal)
     blackBox.ideal = equationsIdeal;      

     -- maybe blackBox.addProperty( equations, gens equationsIdeal )
     blackBox.equations =  gens  equationsIdeal; 


    -- setValuesAt generates 'isZeroAt' and 'jacobianAt', too !  
    blackBox.registerPointProperty("valuesAt",(point)->  (    return  gens sub( equationsIdeal , point);   ));

    -- maybe blackBox.addProperty( jacobian, jacobian gens  equationsIdeal )
    blackBox.jacobian = jacobian gens  equationsIdeal;

    -- needs updatePointProperty, because "jacobianAt" is present. 
    -- Todo: introduce force mode option? (overwriting registered property?)
    -- todo: introduce for all methods a precomposition with checkInputPoint
    blackBox.updatePointProperty( "jacobianAt",

       (point)->
       (   
          jacobianM2MatrixAt := sub( blackBox.jacobian , point);
          return jacobianM2MatrixAt;
       )
    );   

   
     blackBox.setThis(blackBox);
     blackBox.clearInternal();

     return new HashTable from blackBox;
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
   

   point = matrix {{3}};
   try {    IFPBlackBox.jacobianAt(point) } then 
   (
     error("testblackBoxIdeal: jacobianAt should fail due the coefficient ring of the point matrix does not match the ideal coefficient ring and also the ideal coefficients are not integers ");
   )   else
   ();
   point = sub( point, coeffRng ) ;
   try {    IFPBlackBox.jacobianAt(point) } then 
   (    
   ) else
   (
        error("testblackBoxIdeal: jacobianAt should succeed  due the coefficient ring of the point matrix matches the ideal coefficient ring ");
   );  
    IFPBlackBox.ring;
    IFPBlackBox.valuesAt(point) ;
    assert(   IFPBlackBox.isZeroAt( point ) );
   assert( IFPBlackBox.jacobianAt(point)==sub( jacobian IFP,point) );
   assert( IFPBlackBox.valuesAt(point)== gens sub(  IFP, point ) );
)



blackBoxIdealFromEvaluation = method();

blackBoxIdealFromEvaluation(ZZ, Ring, Function) := HashTable => ( numVariables, coeffRing, valuesAt )  ->
(
    blackBox := basicBlackBox(numVariables, coeffRing);

    blackBox.setThis(blackBox);

    --blackBox.setValuesAt ( valuesAt ); --sets isZeroAt, jacobianAt and imageRank

    blackBox.registerPointProperty ("valuesAt", valuesAt ); --sets isZeroAt, jacobianAt and imageRank

    check := ()->
    (
         numVariables :=  blackBox.numVariables;

         point := matrix { apply(numVariables, i-> 0_(blackBox.coefficientRing) ) };
         blackBox.valuesAt( point );
         blackBox.isZeroAt( point );
    );

     check();  
   
    blackBox.setThis(blackBox);-- not necessary?
    blackBox.clearInternal();

    return new HashTable from blackBox ;
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


--blackBoxIdealFromProperties = method();

--blackBoxIdealFromProperties(ZZ, Ring, Function) := HashTable => ( numVariables, coeffRing, pPropertiesAt )  ->
--(
--    blackBox := basicBlackBox( numVariables, coeffRing );
--    blackBox.setThis(blackBox);
--
--    blackBox.setIsZeroAt( (point)->(true) );
--    propertiesAt  := (point)->null;
--
--    --blackBox.propertiesAt = (point)->propertiesAt(point);
   
--    --blackBox.setPropertiesAt ( pPropertiesAt );


--    blackBox.setThis(blackBox);
--    blackBox.clearInternal();

--    return new HashTable from blackBox;
--)

--beginDocumentation()
    
    

doc ///
   Key
        blackBoxIdeal
   Headline
        create a blackbox describing an ideal from an ideal with integer coefficient ring
   Usage   
        blackBoxIdeal(IdealInZZ)
        blackBoxIdealFromEvaluation(rng,evaluationMethod)
        blackBoxIdealFromEvaluation(variableNumber,coeffRing, evaluationMethod)
   Inputs  
        IdealInZZ:Ideal
             with integer coefficients
   Outputs
        : HashTable
             an ideal black box with methods  \break
            \,\, \bullet \, { \tt knownPointProperties() }. \break
            \,\, \bullet \, { \tt registerPointProperty(pointPropertyMethod) }. \break
            where   { \tt pointPropertyMethod } expects to take two parameters: (blackBox,point)  \break
            \,\, \bullet \, { \tt valuesAt(point)},\break
            \,\, \bullet \, { \tt jacobianAt(point) }. \break

            optional:
            \,\, \bullet \, { \tt equations}, \break
            \,\, \bullet \, { \tt unknowns},\break
            \,\, \bullet \, { \tt jacobian}, \break
   Description
        Text
            Creates a blackbox describing an ideal from an ideal with integer coefficient ring
        Example          
        Text
           \break  Example:  create an ideal black box object from an ideal:
        Example          
            RQ := QQ[x];
            IFQ := ideal { 1/3*x^2+1, 1/5*x+2 };        
            IFZ := clearCoeffDenominators(IFQ)
            IFZBlackBox := blackBoxIdeal(IFZ);
        Text
            \break Now access some ideal propeties via the black box interface:
        Example          
            keys IFZBlackBox
            IFZBlackBox.unknowns
            IFZBlackBox.equations
            IFZBlackBox.jacobian
            point := matrix { {1} };
            IFZBlackBox.valuesAt(point)
            IFZBlackBox.jacobianAt(point)            
   Caveat
        does not check if the ideal ring is a quotient ring (not supported)
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
assert (B.jacobianAt(line) == sub(jacobian I,line))

assert (B2.bareJacobianAt(line) == B.bareJacobianAt(line))

///


TEST ///

    bbRankM = blackBoxIdeal( 5 ,ZZ )
    assert(bbRankM.numVariables==5);
    assert(bbRankM.coefficientRing===ZZ);

    assert( bbRankM.imageRank() === null)

    rankMat := (point)->5 

   assert(not  bbRankM.hasPointProperty("rankMat") );

    bbRankM.registerPointProperty("rankMat",rankMat)

    assert(  bbRankM.hasPointProperty("rankMat") );

    point  = matrix {{1,2,3,4,5}};

    assert( rankMat(point) == (bbRankM.pointProperty("rankMat"))(point) );

    rankMatNew := (point)->3
    bbRankM.updatePointProperty("rankMat",rankMatNew)
    assert( rankMatNew(point) == (bbRankM.pointProperty("rankMat"))(point) );
    assert( rankMatNew(point) == (bbRankM.pointProperty(getGlobalSymbol "rankMat"))(point) );


    bbRankMNew = bbRankM.rebuild()
    assert bbRankMNew#?(global rankMat);
    assert bbRankMNew#?("rankMat");

    assert( rankMatNew(point) == bbRankMNew.rankMat(point) );

    assert(bbRankMNew.coefficientRing===ZZ);

    rankMatNew := (point)->4 --also influences bbRankM; because rebuild does not copy; it just exports new registered properties.

    bbRankMNew.updatePointProperty("rankMat",rankMatNew)

    (bbRankM.pointProperty("rankMat"))(point);

    assert( rankMatNew(point) == bbRankMNew.rankMat(point) );
    assert( rankMatNew(point) == bbRankMNew#"rankMat"(point) );
    assert( rankMatNew(point) == (bbRankMNew.pointProperty("rankMat"))(point) );
    assert( rankMatNew(point) == (bbRankMNew.pointProperty(getGlobalSymbol "rankMat"))(point) );
    keys bbRankM
    
    valuesAt := (point)-> matrix {{1,2}};

    bbRankM.registerPointProperty("valuesAt",valuesAt);

    assert(  bbRankM.hasPointProperty("isZeroAt") );

    assert(  bbRankM.hasPointProperty("jacobianAt") );
    assert(  bbRankM.hasPointProperty("bareJacobianAt") );

    assert( bbRankM.imageRank() =!= null)

    assert( bbRankM.imageRank() === 2 )
    bbRankMNew.imageRank()
    
    illegalPoint := matrix {{1,2,3,4,5,6}}; 
   
    try ( bbRankMNew.rankMat(illegalPoint) ) then ( assert(false) ) else ();


    bbRankM = blackBoxIdeal( 5 ,ZZ/7 )

    valuesAt := (point)-> matrix {{1,2}};

    bbRankM.registerPointProperty("valuesAt",valuesAt);

    bbRankM = bbRankM.rebuild();

    point  = sub(point,ZZ/7); 

    bbRankM.valuesAt(point)

    illegalPoint := sub(point,ZZ/2); 

    try ( bbRankM.valuesAt(illegalPoint) ) then ( assert(false) ) else ();

   
 ///

end
-- need test for randomIterator (for fixed error : to less trials  )



