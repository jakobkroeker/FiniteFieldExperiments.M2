



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
       print ("orig blackBox.ring  called");
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

     blackBox.sourceRank = ()->
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
       --point := matrix { apply(blackBox.sourceRank(), i-> 0_(blackBox.coefficientRing()) ) };
       --return #blackBox.valuesAt( point );
   );


   blackBox.setImageRank = (pImageRank)->
   (
      imageRank = pImageRank;
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

     transposedJacobian := jacobian gens equationsIdeal;
     jacobianMatrix := transpose transposedJacobian;

     --transposedJacobian := symbol transposedJacobian;

     -- should be different if the ideal is over Fp

      blackBox.setRing(ring equationsIdeal);
      remove( blackBox, getSymbol "setRing" );

  

     remove( blackBox, getSymbol "setImageRank" );


     blackBox.transposedJacobian=() ->
     (
       return transposedJacobian;
     );

     --blackBox.jacobianMatrix = () ->
     --(
     --  return jacobianMatrix;
     --);

   blackBox.jacobian = () ->
     (
       return jacobianMatrix;
     );

     blackBox.getEquations = ()->
     (
         return gens  equationsIdeal;   
     );

     valuesAt := (point)->
     (
         return transpose gens sub( equationsIdeal , point);   
     );

     blackBox.setValuesAt(valuesAt);
     remove( blackBox, getSymbol "setValuesAt" );
     blackBox.internalRegisterPointProperty( "valuesAt" , blackBox.valuesAt );

     jacobianAt := (point)->
     (
        jacobianMatrixAt:= sub( jacobianMatrix , point);
        --get rid of map degree information
	jacobianMatrixAt = sub( jacobianMatrixAt, ring point);
        jacobianMatrixAt = sub(sub(jacobianMatrixAt, ZZ), ring point);
        return jacobianMatrixAt;
     );

       
    blackBox.jacobianAt= (point)->
    (
        blackBox.checkCoeffRing(point);
        return jacobianAt(point);
    );      
     
   --blackBox.setImageRank(rank image gens equationsIdeal);
    tmppoint := matrix { apply(blackBox.sourceRank(), i-> 0_(blackBox.coefficientRing()) ) };
     blackBox.setImageRank(numRows blackBox.valuesAt( tmppoint ));


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
   assert( IFPBlackBox.jacobian()==transpose jacobian IFP);
   

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
   assert( IFPBlackBox.jacobianAt(point)==sub(transpose jacobian IFP,point) );
   assert( IFPBlackBox.valuesAt(point)==transpose gens sub(  IFP, point ) );
)

blackBoxIdealFromEvaluation  :=method();

blackBoxIdealFromEvaluation( Ring, Function )  := HashTable=> ( rng, valuesAt )->
(

     blackBox := new MutableHashTable from createBasicBlackBox();
 
      blackBox.setRing(rng);
      remove( blackBox, getSymbol "setRing" );

     blackBox.ring = ()->
     (         
        print ("new blackBox.ring  called");
        return rng;
     );

    
     blackBox.setValuesAt(valuesAt);
     remove( blackBox, getSymbol "setValuesAt" );

     blackBox.valuesAt = (point)->
     (
         print ("new blackBox.valuesAt never called");
         return valuesAt( point);   
     );
 
  
     tmppoint := matrix { apply(blackBox.sourceRank(), i-> 0_(blackBox.coefficientRing()) ) };
     blackBox.setImageRank(numRows blackBox.valuesAt( tmppoint ));
     remove( blackBox, getSymbol "setImageRank" );

     blackBox.isZeroAt = (point)->
     (
         return blackBox.valuesAt(point)==0;   
     );

     check := ()->
     (
         sourceRank :=  blackBox.sourceRank();

         point := matrix { apply(sourceRank, i-> 0_(blackBox.coefficientRing()) ) };
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

        jacobianMatrixAt := mutableMatrix( blackBox.coefficientRing() ,  blackBox.imageRank(  ), #unknowns );
        for unknownIdx  in 0..(#unknowns-1) do
        (
             newpoint := new MutableMatrix from sub( point,epsRng );
   
             newpoint_(unknownIdx,0) = newpoint_(unknownIdx,0)+eps;
             valueVec := blackBox.valuesAt( matrix newpoint );  
             for rowIdx in 0..numRows valueVec-1 do
             (
                coordinateValue := last coefficients (valueVec_(rowIdx,0), Monomials=>{1 , eps } );
                if ( not (coordinateValue)_(0,0) ==0) then error("error in jacobianAt. please contact the developers");
                jacobianMatrixAt_(rowIdx,unknownIdx) = sub( (coordinateValue)_(1,0) , blackBox.coefficientRing())  ;
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
    a = getSymbol a;
    
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

    assert( evaluation.sourceRank() ==evalBlackBox.sourceRank() );

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
debug padicLift
padicLiftProtect()
testBlackBoxIdeal()
///

TEST ///
debug padicLift
padicLiftProtect()
testBlackBoxIdealFromEvaluation()
///

end
---

restart
--load"varietyBlackBoxes.m2"
loadPackage"padicLift"
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
	  r = random(K^1,K^4)
	  assert (B2.isZero(r) == B.isZero(r))
	  assert (B2.valuesAt(r) == B.valuesAt(r))
	  assert (B2.jacobianAt(r) == B.jacobianAt(r))
     ))

assert B2.isZeroAt(line)


restart
loadPackage"padicLift"
