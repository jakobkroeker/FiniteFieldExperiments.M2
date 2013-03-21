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

export{
    clearCoeffDenominators,
    blackBoxIdeal,
    blackBoxIdealFromEvaluation,
    blackBoxIdealFromProperties
}


idealBlackBoxesProtect = ()->
(
   
    protect jacobianAt;
    protect jacobianMatrix;
    protect transposedJacobianAt;
    protect transposedJacobian;
    protect valuesAt;
    protect unknownIsValid;
    protect getUnknowns;
    protect createBlackBoxIdeal;
    protect getEquations;
    protect numVariables;
    protect imageRank;
    protect isZeroAt;
    protect getPointProperties;
    protect setRing;
    protect registerPointProperty;
    protect internalRegisterPointProperty;
    protect setValuesAt;
    protect setImageRank;
    protect checkCoeffRing;
    protect deduceImageRank;
    protect properties;

)

--todo: fix dublicate code,  -  padicLiftProtect and padicLiftExport
idealBlackBoxesExport = ()->
(
    exportMutable( jacobianAt);
    exportMutable( jacobianMatrix);
    exportMutable( transposedJacobian);
    exportMutable( transposedJacobianAt);
    exportMutable( valuesAt);
    exportMutable( unknownIsValid);    
    exportMutable( getUnknowns);
    exportMutable( createBlackBoxIdeal);
    exportMutable( getEquations);                
    exportMutable( numVariables);  
    exportMutable( imageRank);
    exportMutable( isZeroAt);      
    exportMutable( getPointProperties);  
    exportMutable( setRing);  
    exportMutable( registerPointProperty);  
    exportMutable( internalRegisterPointProperty);
    exportMutable( setValuesAt);    
    exportMutable( setImageRank );
    exportMutable(  checkCoeffRing);
    exportMutable( deduceImageRank );
    exportMutable( properties );
)




-- swith between protect and export - both are not possible!

--idealBlackBoxesProtect() -- protect the symbols for checking package correctness: no symbol variables should be overwritten by accident!
idealBlackBoxesExport(); -- export the symbols to make the package work 


beginDocumentation()


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
clearCoeffDenominators = method();

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




createBasicBlackBox = () ->
(
    
   blackBox := new MutableHashTable;
      

   pointProperties := new HashTable  ;

   blackBox.getPointProperties=()->
   (  
      return pointProperties;
   );

    
     rng := null;

     blackBox.ring = () ->
     (
        --error "error:  blackBox.ring() not implemented";
        return rng;
     );

  
     blackBox.setRing = (parrng) ->
     ( 
        rng =  parrng;
     );

     blackBox.coefficientRing = () ->
     (
       return  coefficientRing blackBox.ring();
     );

     blackBox.checkCoeffRing = (point)->
     (
        errorMsg := "  ideal is defined over "| toString  blackBox.coefficientRing() | "and not over " | toString ring point |"!" ;
        if blackBox.coefficientRing() =!= ZZ  
           and ring point =!= blackBox.coefficientRing()
           then 
	   (
		try ( if coefficientRing ring point=!= blackBox.coefficientRing() then error errorMsg )  then () else ( error (errorMsg) );
           );
     );

     blackBox.getUnknowns =()->
     (
         return gens blackBox.ring();
     );

     blackBox.unknownIsValid = (unknown)->
     (
        if not ( blackBox.ring() === ring unknown) then 
        ( 
            print "the unknown is not element of the equations ideal ring";
	    return false;
        );
        return true;
     );

     blackBox.numVariables = ()->
     (
        --return rank source  blackBox.getUnknowns();
        return #(blackBox.getUnknowns());
     );

    blackBox.registerPointProperty = method();
   blackBox.registerPointProperty(String, Function, Function) := Thing => ( propertyName, propertyMethod, preconditionsTest )->
   (

      preconditionsTest(blackBox); -- e.g. check for method presence or test an example.
      
      pointProperties = new MutableHashTable from pointProperties;

      propertySymbol := getSymbol propertyName;
      if  ( not  pointProperties#?propertySymbol  and  not  pointProperties#?propertyName 
      and not  blackBox#?propertySymbol          and   not  blackBox#?propertyName  ) then 
      (
         pointProperties#propertySymbol = (point )-> (return  propertyMethod( blackBox,point); );
         pointProperties#propertyName   = (point )-> (return  propertyMethod( blackBox,point); );
         blackBox#propertySymbol        = (point )-> (return  propertyMethod( blackBox,point); );
         blackBox#propertyName          = (point )-> (return  propertyMethod( blackBox,point); );
         return new HashTable from blackBox;
      )
      else "error: method "| propertyName |" is already registered";
 
      pointProperties = new HashTable from pointProperties;

     -- error "not yet implemented"; -- idea :  propertyMethod takes two parameters; first one is the ideal black box, second one is point
   );

   -- how to hide this method?
   -- propertyMethod is expected to have two parameters: (blackBox,point)
   blackBox.registerPointProperty(String, Function) := Thing =>( propertyName, propertyMethod )->
   (
      blackBox.registerPointProperty( propertyName, propertyMethod, ()->() )
   );

   -- how to hide this method?
   blackBox.internalRegisterPointProperty = method();
   blackBox.internalRegisterPointProperty (String, Function, Function) := Thing => ( propertyName, propertyMethod, preconditionsTest )->
   (

      preconditionsTest(blackBox); -- e.g. check for method presence or test an example.

      pointProperties = new MutableHashTable from pointProperties;
      propertySymbol := getSymbol propertyName;
      if not  pointProperties#?propertySymbol  and  not  pointProperties#?propertyName then 
      (
         pointProperties#propertySymbol = (point )-> (return  propertyMethod(point); );
         pointProperties#propertyName = (point )-> (return  propertyMethod(point); );
      )
      else "error: method "| propertyName |" is already registered";
 
      pointProperties = new HashTable from pointProperties;
      return null;
   );

   blackBox.internalRegisterPointProperty(String, Function) := Thing => ( propertyName, propertyMethod )->
   (
      blackBox.internalRegisterPointProperty( propertyName, propertyMethod, (varietyIdeal)->() );
      return null;
   );

   valuesAt := null;

   blackBox.valuesAt = (point) ->
   (
      blackBox.checkCoeffRing(point);
      return valuesAt(point);
   );

   imageRank := null;

   blackBox.imageRank = ()->
   (
       return imageRank;
       --point := matrix { apply(blackBox.numVariables(), i-> 0_(blackBox.coefficientRing()) ) };
       --return #blackBox.valuesAt( point );
   );


   blackBox.setImageRank = (pImageRank)->
   (
      imageRank = pImageRank;
   );

     blackBox.deduceImageRank = ()->
     (
        computed  := false; 
        maxTrials := 100;
        currTrial := 0;
        rng := blackBox.coefficientRing();
        imageRank := null;
        while imageRank===null and currTrial<maxTrials do
        (
          try (
              tmppoint := matrix random(rng^1,rng^(blackBox.numVariables()) );
              valuesMatrix := blackBox.valuesAt( tmppoint );
              print valuesMatrix;
              assert (numRows valuesMatrix==1);
              imageRank = numColumns valuesMatrix;

          )         then ( computed=true ) else();
               
           
          currTrial = currTrial+1;
        );
        if not computed then error "error: failed to deduce image dim";
        return imageRank;
     );


    blackBox.setValuesAt = (pValuesAt) ->
     ( 
        valuesAt =  pValuesAt;
     );


   blackBox.isZeroAt = (point)->
   (
       return blackBox.valuesAt(point)==0;   
   );


   blackBox.internalRegisterPointProperty( "isZeroAt" , blackBox.isZeroAt );

   return new HashTable from blackBox;
)






-- todo: how to check, if 'ring equationsIdeal' is not a quotient ring?
blackBoxIdeal  = (equationsIdeal)->
(
     blackBox := new MutableHashTable from createBasicBlackBox();

     jacobianM2Matrix := jacobian gens equationsIdeal;


     -- should be different if the ideal is over Fp

      blackBox.setRing(ring equationsIdeal);
      remove( blackBox, getSymbol "setRing" );


     blackBox.jacobian = () ->
     (
       return jacobianM2Matrix;
     );

     blackBox.getEquations = ()->
     (
         return gens  equationsIdeal;   
     );

     valuesAt := (point)->
     (
         return gens sub( equationsIdeal , point);   
     );

     blackBox.setValuesAt(valuesAt);
     remove( blackBox, getSymbol "setValuesAt" );
     blackBox.internalRegisterPointProperty( "valuesAt" , blackBox.valuesAt );

     jacobianAt := (point)->
     (
        jacobianM2MatrixAt:= sub( jacobianM2Matrix , point);
        --get rid of map degree information
	    jacobianM2MatrixAt = sub( jacobianM2MatrixAt, ring point);
        jacobianM2MatrixAt = sub(sub(jacobianM2MatrixAt, ZZ), ring point);
        return jacobianM2MatrixAt;
     );

       
    blackBox.jacobianAt= (point)->
    (
        blackBox.checkCoeffRing(point);
        return jacobianAt(point);
    );      
     
     blackBox.setImageRank( blackBox.deduceImageRank() );
     remove( blackBox, getSymbol "setImageRank" );


     blackBox.internalRegisterPointProperty( "jacobianAt" , blackBox.jacobianAt );
  
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
   assert( IFPBlackBox.getUnknowns()=={x} );
   assert( IFPBlackBox.getEquations()==gens IFP);
   assert( IFPBlackBox.jacobian()== jacobian IFP);
   

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
    IFPBlackBox.ring();
    IFPBlackBox.valuesAt(point) ;
    assert(   IFPBlackBox.isZeroAt( point ) );
   assert( IFPBlackBox.jacobianAt(point)==sub( jacobian IFP,point) );
   assert( IFPBlackBox.valuesAt(point)== gens sub(  IFP, point ) );
)

blackBoxIdealFromEvaluation = method();

blackBoxIdealFromEvaluation( Ring, Function )  := HashTable=> ( rng, valuesAt )->
(

     blackBox := new MutableHashTable from createBasicBlackBox();
 
      blackBox.setRing(rng);
      remove( blackBox, getSymbol "setRing" );

     blackBox.ring = ()->
     (         
        return rng;
     );

    
     blackBox.setValuesAt(valuesAt);
     remove( blackBox, getSymbol "setValuesAt" );

     blackBox.valuesAt = (point)->
     (
         return valuesAt( point);   
     );
 
    

     blackBox.setImageRank( blackBox.deduceImageRank() );
     remove( blackBox, getSymbol "setImageRank" );

     blackBox.isZeroAt = (point)->
     (
         return blackBox.valuesAt(point)==0;   
     );

     check := ()->
     (
         numVariables :=  blackBox.numVariables();

         point := matrix { apply(numVariables, i-> 0_(blackBox.coefficientRing()) ) };
         blackBox.valuesAt( point );
         blackBox.isZeroAt( point );
     );

     check();

     -- here comes the magic:
     blackBox.jacobianAt = (point)->
     (
       
        unknowns := blackBox.getUnknowns() ;      
        
        if (not ( blackBox.valuesAt( point )==0))  then  error("point does not belong to the ideal ! ");


        eps := null;
        eps = symbol eps;
        epsRng := ( blackBox.coefficientRing() )[eps]/eps^2;
        eps = (gens epsRng)#0;

        jacobianMatrixAt := mutableMatrix( blackBox.coefficientRing(),  #unknowns, blackBox.imageRank(  ) );
        for unknownIdx  in 0..(#unknowns-1) do
        (
             newpoint := new MutableMatrix from sub( point,epsRng );
   
             newpoint_(unknownIdx,0) = newpoint_(unknownIdx,0)+eps;
             valueVec := blackBox.valuesAt( matrix newpoint );  
             for equationIdx in 0..numColumns valueVec-1 do
             (
                coordinateValue := last coefficients (valueVec_(0,equationIdx), Monomials=>{1 , eps } );
                if ( not (coordinateValue)_(0,0) ==0) then error("error in jacobianAt. please contact the developers");
                jacobianMatrixAt_(unknownIdx,equationIdx) = sub( (coordinateValue)_(1,0) , blackBox.coefficientRing())  ;
             );
        );
        return matrix jacobianMatrixAt;
     );

   
    blackBox.internalRegisterPointProperty("jacobianAt" , blackBox.jacobianAt );
    blackBox.internalRegisterPointProperty("valuesAt" , blackBox.valuesAt );
  

    return  new HashTable from blackBox;
)

blackBoxIdealFromEvaluation(ZZ, Ring, Function) := HashTable => ( numVariables, coeffRing, valuesAt )  ->
(
    assert ( numVariables>0 );
    a := null;
    a = symbol a;
    -- todo: add getNumVariables;

    rng := coeffRing[a_1..a_numVariables];
    return blackBoxIdealFromEvaluation(rng,valuesAt);
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
    --evaluation := createBasicBlackBox();
  
    evalBlackBox := blackBoxIdealFromEvaluation ( evaluation.ring(), evaluation.valuesAt );

    point := matrix {{3_(ZZ/7)}} ;
    assert( evaluation.isZeroAt( point ) );
    assert( evaluation.getUnknowns()=={x} );
    assert( evaluation.valuesAt( point ) == evalBlackBox.valuesAt( point ) );
    assert( evaluation.jacobianAt( point ) == evalBlackBox.jacobianAt( point ) );

    assert( evaluation.numVariables() ==evalBlackBox.numVariables() );

    outerPoint := matrix {{2_(ZZ/7)}} ;

    assert( evaluation.valuesAt( outerPoint ) == evalBlackBox.valuesAt( outerPoint ) );

    assert( not evaluation.isZeroAt( outerPoint ) );

    assert( evaluation.coefficientRing() === coeffRng);

    assert( evaluation.unknownIsValid ( ( gens ring IFP )#0 ));

    y  := null;    y  = symbol y;
    rngy := ZZ/7[y];
    y = (gens rng)#0;

    assert( not evaluation.unknownIsValid (  y ) );

);


blackBoxIdealFromProperies = method();

blackBoxIdealFromProperies(ZZ, Ring, Function) := HashTable => ( numVariables, coeffRing, properties )  ->
(
    assert ( numVariables>0 );
    B := new MutableHashTable;
    B.numVariables = ()->numVariables;
    B.coefficientRing = ()->coeffRing;
    B.properties = properties;
    return B
)

    
    


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
            \,\, \bullet \, { \tt getPointProperties() }. \break
            \,\, \bullet \, { \tt registerPointProperty(pointPropertyMethod) }. \break
            where   { \tt pointPropertyMethod } expects to take two parameters: (blackBox,point)  \break
            \,\, \bullet \, { \tt valuesAt(point)},\break
            \,\, \bullet \, { \tt jacobianAt(point) }. \break

            optional:
            \,\, \bullet \, { \tt getEquations()}, \break
            \,\, \bullet \, { \tt getUnknowns()},\break
            \,\, \bullet \, { \tt jacobian()}, \break
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
            IFZBlackBox.getUnknowns()
            IFZBlackBox.getEquations()
            IFZBlackBox.jacobian()
            point := matrix { {1} };
            IFZBlackBox.valuesAt(point)
            IFZBlackBox.jacobianAt(point)            
   Caveat
        does not check if the ideal ring is a quotient ring (not supported)
///


TEST ///
debug idealBlackBoxes
idealBlackBoxesProtect()
testClearCoeffDenominators()
///



TEST ///
debug idealBlackBoxes
idealBlackBoxesProtect()
testNestedRingCoeffsLCMDenominator()
///
         

TEST ///
debug idealBlackBoxes
idealBlackBoxesProtect()
testTensoredClearCoeffDenominators()
///

TEST ///
debug idealBlackBoxes
idealBlackBoxesProtect()
testBlackBoxIdeal()
///

TEST ///
debug idealBlackBoxes
idealBlackBoxesProtect()
testBlackBoxIdealFromEvaluation()
///

end
---

restart
--load "idealBlackBoxes.m2"
loadPackage"idealBlackBoxes"
apropos "blackBox"

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
rank(B.jacobianAt smoothPoint) == codim(I,Generic=>true)
rank(B.jacobianAt singularPoint) < codim(I,Generic=>true)
B.isZeroAt(smoothPoint)
B.isZeroAt(singularPoint)
B.isZeroAt(line)
not B.isZeroAt(offPoint)

prime = 11
K = ZZ/prime
B.isZeroAt(sub(smoothPoint,K))
B.isZeroAt(sub(offPoint,K))

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

B2 = blackBoxIdealFromEvaluation(4,ZZ,evalLinePlusConic)

apply(100,i->(
	  r = random(K^1,K^4);
	  assert (B2.isZeroAt(r) == B.isZeroAt(r));
	  assert (B2.valuesAt(r) == B.valuesAt(r));
      	  if B2.isZeroAt(r) then 
	      assert (B2.jacobianAt(r) == B.jacobianAt(r));
     ))

assert B2.isZeroAt(line)
B2.jacobianAt(line)
B.jacobianAt(line)

-- sollte gehn
assert (B.jacobianAt(line) == sub(jacobian I,line))
-- eventuell schwieriger
assert (B2.jacobianAt(line) == B.jacobianAt(line))



line


restart
loadPackage"idealBlackBoxes"
