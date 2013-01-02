

adjustRootPairingTolerance=  ( tolerance, rootList ) ->
(
   localTolerance := tolerance;
   numRoots := #rootList;
   bitPrecision := precision rootList#0;
   denom := value ("3.0" |"p" | toString bitPrecision);

   for row in 0..numRoots-1 do
   (
     for col in (row+1)..numRoots-1 do
     (
       toleranceCandidate:= abs( rootList#row - rootList#col / denom );
       if toleranceCandidate <abs(localTolerance) then
           localTolerance  =toleranceCandidate;
     );
   );
   return localTolerance;
)


rootCompatibilityMatrixRowsValid= ( compatibiltyMatrix, exact )->
(   
     rowNonzeroCount := apply(numRows compatibiltyMatrix, rowNum-> sum(apply(flatten entries compatibiltyMatrix^{rowNum} ,i-> if i>0 then 1 else 0)) );
 
     for nonzeroCount in rowNonzeroCount do
     (
        if nonzeroCount>1 and exact then
            return false;
 
        if nonzeroCount<1 then
          return false;
    );
    return true;
);

testRootCompatibilityMatrixRowsValid=()->
(
	compatibiltyMatrix := matrix {{ 0,1},{2,0},{0,3}};
	assert( rootCompatibilityMatrixRowsValid(compatibiltyMatrix,false) );
	assert( rootCompatibilityMatrixRowsValid(compatibiltyMatrix,true) );

	compatibiltyMatrix = matrix{{2,1},{2,0},{0,3}};
	assert( rootCompatibilityMatrixRowsValid(compatibiltyMatrix,false) );

	compatibiltyMatrix = matrix{{2,1},{2,0},{0,3}};
	assert(  not rootCompatibilityMatrixRowsValid(compatibiltyMatrix,true) );

	compatibiltyMatrix = matrix {{2,1},{0,0},{0,3}};
	assert(  not rootCompatibilityMatrixRowsValid(compatibiltyMatrix,true) );
	assert(  not rootCompatibilityMatrixRowsValid(compatibiltyMatrix,false) );
)
 

TEST ///
debug padicLift
padicLiftProtect()
  testRootCompatibilityMatrixRowsValid();
///


isValidRootCompatibility = ( Mat, combinedRootsCount )->
(
 assert(char ring Mat == 0 );
 
    mathchedRoots := set flatten entries Mat  ;
    mathchedRoots = mathchedRoots - {0} ;
 
    -- for each combined root there should be a existing compatibility:
    if #mathchedRoots!=combinedRootsCount then
    (
         --Info(InfoHMACRootPairing, 1, "---------root compatibility warning: Size(mathchedRoots)combinedRootsCount, problem with error tolerance?" );
        return false ;
    );

     if not rootCompatibilityMatrixRowsValid ( Mat, false) or
        not rootCompatibilityMatrixRowsValid ( transpose Mat,false) then
     (
               --Info(InfoHMACRootPairing, 1,"-------------root compatibility warning: compatibility not given; problem with error tolerance ?");
        return false ;
     );
    return true;
)


isValidRootCompatibilityTest=()->
(
    Mat := matrix  {{1,2},{1,4},{5,6 }};
    assert(not isValidRootCompatibility (Mat,6) );
    Mat  =  matrix {{1,2},{3,4},{5,6}};
    assert( isValidRootCompatibility (Mat,6) );
    assert( isValidRootCompatibility (Mat,6) );
    
     Mat =  matrix {{1,0},{3,0},{2,0}};
     assert( not isValidRootCompatibility (Mat,3) );
     
      Mat =  matrix {{1,0},{0,3},{2,0}};
     assert( isValidRootCompatibility (Mat,3) );
);

TEST ///
debug padicLift
padicLiftProtect()
  isValidRootCompatibilityTest();
///


computeWeakRootCompatibility=( firstPolRoots, secondPolRoots, combinedPolRoots, operation, maxTolerance)->
(
    localTolerance := maxTolerance;
    
    assert(#firstPolRoots >= #secondPolRoots  );

    if not #firstPolRoots >=#secondPolRoots  or
       not #combinedPolRoots == max( #firstPolRoots ,  #secondPolRoots ) then
          return null;

  
    numRoots := #firstPolRoots ;
    compatibiltyMatrix :=  mutableMatrix (ZZ,numRoots,#secondPolRoots)  ;

    extendedCompatibilityMatrix :=  mutableMatrix (ZZ,numRoots,#secondPolRoots)  ;
    
    --# tolerance is after adjusting smaller or equal to the minimal distance between two roots for each root list
    localTolerance = adjustRootPairingTolerance ( localTolerance, firstPolRoots );
    localTolerance = adjustRootPairingTolerance ( localTolerance, secondPolRoots );
    localTolerance = adjustRootPairingTolerance ( localTolerance, combinedPolRoots );


    if abs( localTolerance)==0 then
    ( 
         --Info(InfoHMACRootPairing, 1, "ComputeWeakRootCompatibility@HMAC: pairing tolerance is zero ");
         return null;
    );
   
    
	for row in 0..numRoots-1 do
	for col in 0..#(secondPolRoots)-1 do
	for i in 0..numRoots-1 do
		if abs( operation (firstPolRoots#row , secondPolRoots#col )- combinedPolRoots#i ) <abs(localTolerance) then
		(
		    compatibiltyMatrix_(row,col) = 1;
		    extendedCompatibilityMatrix_(row,col) = i+1;
		);
    
    if not (rank matrix compatibiltyMatrix == #secondPolRoots ) then
        return null;

    rowSums := apply( numRows compatibiltyMatrix, i-> sum(flatten entries (matrix compatibiltyMatrix)^{i} ) );
    for entry in rowSums do
    (
        if not entry==1 then
          return null;
    );

    return matrix compatibiltyMatrix;
)

computeRootCompatibility=
 ( firstPolRoots, secondPolRoots, combinedPolRoots, operation, maxTolerance)->
(   
    localTolerance := maxTolerance;
    
    if not #combinedPolRoots >= max( # firstPolRoots , #secondPolRoots  ) then
            --Info(InfoHMACRootPairing, 1, "ComputeRootCompatibility: Error: Size(combinedPolRoots)<Maximum( Size(firstPolRoots), Size(secondPolRoots) )" );
          return null;
 

    numRoots := #firstPolRoots ;
    compatibilityMatrix :=  mutableMatrix (ZZ,numRoots,#secondPolRoots)  ;

    simpleCompatibiltyMatrix :=  mutableMatrix (ZZ,numRoots,#secondPolRoots)  ;
  
    
    localTolerance = adjustRootPairingTolerance ( localTolerance, firstPolRoots );
    localTolerance = adjustRootPairingTolerance ( localTolerance, secondPolRoots );
    localTolerance = adjustRootPairingTolerance ( localTolerance, combinedPolRoots );


    if abs( localTolerance)==0 then
    ( 
         --Info(InfoHMACRootPairing, 1, "ComputeWeakRootCompatibility@HMAC: pairing tolerance is zero ");
         return null;
    );
   
    combinedRootsMatched :=  new MutableList;
    for row in   0..#firstPolRoots-1  do
    for col in   0..#secondPolRoots-1  do
    for i   in   0..#combinedPolRoots-1 do
    (
        if abs( operation (firstPolRoots#row, secondPolRoots#col )- combinedPolRoots#i ) < abs(localTolerance) then
        (
            combinedRootsMatched#i = 1;
            compatibilityMatrix_(row,col) = i+1;
            simpleCompatibiltyMatrix_(row,col) = 1;
        );
    );

    compatibilityMatrix = matrix compatibilityMatrix;
    -- Info(InfoHMACRootPairing, 2, "compatibilityMatrix");
    -- Info(InfoHMACRootPairing, 2, String(compatibilityMatrix) );

    if not isValidRootCompatibility ( compatibilityMatrix, #combinedPolRoots  ) then
    (
         --Info(InfoHMACRootPairing, 1, "--------------ComputeRootCompatibility: probably a problem with error tolerance....." );
        return null;
    );

    return compatibilityMatrix;
)


testComputeRootCompatibility=()->
(
  
    changeFloatPrecision := (floatList  )->
    (
          desiredPrecision := 1000;
	  return apply(floatList, entry-> toRR(desiredPrecision, entry)); 
    );

    firstPolRoots := changeFloatPrecision( { 0.03, 34.0, 10.0 } );

    secondPolRoots   := changeFloatPrecision( { 5.03, 4.0, 1.0 } );
    combinedPolRoots := changeFloatPrecision( { 4.03, 11.0, 39.02 } );
    
    operation := (a,b) ->(return a+b; );
    
    maxTolerance :=  0.001 ;
    compatibilityMatrix := computeWeakRootCompatibility  ( firstPolRoots, secondPolRoots, combinedPolRoots, operation, maxTolerance );
    assert(  compatibilityMatrix == null);

 
     maxTolerance=  0.02 ;
    compatibilityMatrix  = computeWeakRootCompatibility ( firstPolRoots, secondPolRoots, combinedPolRoots, operation,  maxTolerance );
    assert(  compatibilityMatrix== matrix {{ 0, 1, 0},{1, 0, 0},{0, 0, 1 }} );

    firstPolRoots = changeFloatPrecision( {4.0, 10.0 } );
    secondPolRoots = changeFloatPrecision( { 5.0 } );
    combinedPolRoots = changeFloatPrecision( { 9.0, 15.0} );
  
    maxTolerance =  0.02 ;
    compatibilityMatrix  = computeRootCompatibility ( firstPolRoots, secondPolRoots, combinedPolRoots, operation, maxTolerance  );
    assert(  compatibilityMatrix == matrix {{1},{2}} );
)

TEST ///
debug padicLift
padicLiftProtect()
  testComputeRootCompatibility();
///


-------------------------------------deprecated-------------------------------------------

-- todo: operation optional. 
-- todo: increase tolerance
computeRootCompatibilityOld = method (Options=>{tolerance=>0.0001});
computeRootCompatibilityOld (List,List,List,FunctionClosure) := Matrix=> opts->(firstPolRoots,secondPolRoots, combinedPolRoots,operation)->
(
    localTolerance:=opts.tolerance;
    --assert(#firstPolRoots==#secondPolRoots);
    assert(#firstPolRoots>=#secondPolRoots);
    assert(#combinedPolRoots==max(#firstPolRoots, #secondPolRoots));
    --
    numRoots:=#firstPolRoots;
    compatibiltyMatrix := mutableMatrix(ZZ,numRoots,#secondPolRoots);
    --
    permutationVector  := mutableMatrix(ZZ,numRoots,1);
    --
    for row in 0..numRoots-1 do
    for col in (row+1)..numRoots-1 do
         if (abs(   (firstPolRoots#row-firstPolRoots#col))<localTolerance) then 
            localTolerance=(firstPolRoots#row-firstPolRoots#col);
    --print "tolerance:";
    --print  localTolerance;
    assert( localTolerance!=0 );
    --
    for row in 0..numRoots-1 do
    for col in 0..#secondPolRoots-1 do
    for i in 0..numRoots-1 do
    (
        print  (abs( operation (firstPolRoots#row, secondPolRoots#col)- combinedPolRoots#i));
       -- if ( abs( operation (firstPolRoots#row, secondPolRoots#col)- combinedPolRoots#i ) < opts.tolerance) then
        if ( abs( operation (firstPolRoots#row, secondPolRoots#col)- combinedPolRoots#i ) <localTolerance) then
        (
            compatibiltyMatrix_(row,col) = 1;
            permutationVector_(row,0) = col;
        );
    );
    --
    if ( (rank matrix compatibiltyMatrix) != #secondPolRoots) then error "failed to pair roots...";

    rowSums := apply(entries matrix compatibiltyMatrix,i->sum i);
    apply( rowSums, i->assert(i==1) );
    -- todo: test: in jeder Zeile ein Eintrag!

    return (matrix compatibiltyMatrix);
)


--------------------------------------sandbox----------------------

isPermutation= method();
isPermutation(List)  :=Boolean=>( permutation)->
(
    try    (
        apply(permutation,el->assert(class el===ZZ));
        sortedpermutation:=sort permutation;
        apply(#sortedpermutation,idx->assert(sortedpermutation#idx==idx));    )
    then
        return true
    else
    return false;
)

invert = method();
invert (List)  :=List=>( permutation)->
(
    assert(isPermutation permutation);
    result := new MutableList;
    for pos in 0..#permutation-1 do
        result#(permutation#pos)=pos;
    return new List from result;
)

permute = method();
permute (Thing,List) :=List=>( srcList, permutation )->
(
    assert(isPermutation permutation);
    -- todo: check isPermutation
    --result:=new MutableList;
    result := apply(permutation, pos->srcList#pos );   
    return new (class srcList) from result;
)


