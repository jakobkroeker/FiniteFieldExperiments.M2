restart

--uninstallPackage "BlackBoxIdeals"
--installPackage "BlackBoxIdeals"
loadPackage("BlackBoxIdeals",Reload=>true)

check BlackBoxIdeals;


testRegisterPointProperty=()->
(
   x  := null;
   x  = symbol x;
   rng := ZZ/7[x];
   coeffRng := coefficientRing rng;
   x = (gens rng)#0;

   RP := ZZ/7[x];
   IFP := ideal { 3*x^2+1, 5*x-1 };        
   bb := blackBoxIdeal( IFP );

   pointProperty = (point)->
   (
       return (1,1);
   );
   
   bb.registerPointProperty( "pointProperty", pointProperty ); -- new properties now available via bb.getPointProperties'
   bb.registerPointProperty( "pointProperty2", pointProperty );
   
   -- make properties available via 'bb.$property':
   bb = bb.rebuild(); 

   assert (bb#?"pointProperty")
   assert (bb#?(symbol  pointProperty))
   assert (bb#?(getSymbol  "pointProperty"))

   assert (bb#?"pointProperty2")
   assert (bb#?(symbol  pointProperty2))
   assert (bb#?(getSymbol  "pointProperty2"))

  

)

methods Symbol
dictionary pointProperty

apply( sort keys bb, key->(key,dictionary key))

2
3+3

restart
loadPackage "BlackBoxIdeals"
loadPackage("BlackBoxIdeals",Reload=>true)

    bbRankM = blackBoxIdeal( 5 ,ZZ )
    assert(bbRankM.numVariables==5);
    assert(bbRankM.coefficientRing===ZZ);

    assert( bbRankM.imageRank() === null)

    rankMat := (point)->5 

   assert(not  bbRankM.hasPointProperty("rankMat") );

    bbRankM.registerPointProperty("rankMat",rankMat)

    assert(  bbRankM.hasPointProperty("rankMat") );

    assert (1 == # select( bbRankM.knownPointProperties(), 
                          (key)->(key==="rankMat"))             );
    bbRankM.pointProperty("rankMat")

    point := matrix {{1,2,3,4,5}};

    (bbRankM.pointProperty("rankMat"))(point);


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

    bbRankMNew.registerPointProperty("valuesAt",valuesAt);
    assert( bbRankMNew.imageRank() === 2 )
    illegalPoint := matrix {{1,2,3,4,5,6}}; 
  bbRankMNew.rankMat(illegalPoint)
    try ( bbRankMNew.rankMat(illegalPoint) ) then ( assert(false) ) else {}
     illegalPoint := sub(point,ZZ/2); 
 point = matrix {{1,2,3,4,5}};
     
    bbRankM = blackBoxIdeal( 5 ,ZZ/7 )

    valuesAt := (point)-> matrix {{1,2}};

    bbRankM.registerPointProperty("valuesAt",valuesAt);
    bbRankM = bbRankM.rebuild()
    point = sub(point,ZZ/7); 
   
    bbRankM.valuesAt(point)

    illegalPoint := sub(point,ZZ/2); 

    try ( bbRankM.valuesAt(illegalPoint) ) then ( assert(false) ) else ();

 bbRankM.valuesAt(illegalPoint)
