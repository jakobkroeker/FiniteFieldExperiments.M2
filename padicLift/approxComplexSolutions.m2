


approxComplexSolutionsOld = method (Options=>{"options"=>new LiftOptions,"decimalPrecision"=>10});

approxComplexSolutionsOld (Ideal, Matrix, List) := opts->(inputIdeal, solutionPoint, unknownList)->
(

    if (not pariGpIsPresent() ) then 
       error "please install Pari/GP or build Macaulay2 with pari (add 'pari' to the --enable-build-libraries configure parameter. ) to use 'approxComplexSolutions'";

    --localDecimalPrecision := opts#"options"#"decimalPrecision";
    localDecimalPrecision := opts#"decimalPrecision";

    rootCalculator := opts#"options"#"rootCalculator";

    characteristic:=char ring solutionPoint;

    
    (minimalPolynomialsTable, liftInfo) := computeMinPolys(inputIdeal, solutionPoint, unknownList, 
                                            "options" => opts#"options");

    pairedRootRootList := new MutableList;
    if (liftInfo#"requiredLatticeDimension"=!=null) then 
    (
        rootListList := apply(unknownList, unknown-> rootCalculator(minimalPolynomialsTable#unknown#1, localDecimalPrecision) );   
    
        pairedRootRootList#0 = rootListList#0;
        operationRootList:={};

        operationAdd := (a,b)->(return a+b;);
        operationSub := (a,b)->(return a-b;);
        operationRnd := (a,b)->(return a+(random(ZZ/characteristic))*b;);
        --operationInputList:={operationAdd,operationRnd,operationSub};
        operationInputList := apply(characteristic-1, num-> ((a,b)->(return a+((num+1)_(ring b))*b;)));

        operationUsedList:={ };

        for pos in 1..#unknownList-1 do
        (
            --print "pos";
            --print pos;
            
             -- todo: check if pairing succeeded for each pos 
            for operation in operationInputList do
            (
                unknown :=  operation ((unknownList)#0 , (unknownList)#pos );
                --print "computeMinPolys";
                    (compatibilityMinimalPolynomials, compatLiftInfo ) := computeMinPolys(inputIdeal, solutionPoint, {unknown}, 
                                    "options" => opts#"options");
   
                if (compatLiftInfo#"requiredLatticeDimension"=!=null) then 
                (
                    unknownRoots := computeRootsWithGP( compatibilityMinimalPolynomials#unknown#1, localDecimalPrecision);
                    if (#unknownRoots != # (rootListList#0)) then continue;
                    --
         
                     firstPolRoots :=  flatten rootListList#0;
                     secondPolRoots := flatten rootListList#pos;
                     combinedPolRoots := flatten unknownRoots;

                    compMatrix := computeWeakRootCompatibility( firstPolRoots,secondPolRoots, combinedPolRoots  ,operation, opts#"options"#"maxPairingTolerance" );
                    if compMatrix===null then continue;  -- todo: do not have a test for this!              

                    modSrcRootList := matrix {rootListList#pos};
                    modDstRootList := compMatrix*(transpose modSrcRootList);
                    operationUsedList = operationUsedList | {operation};
                    pairedRootRootList#pos = flatten entries modDstRootList;
                    break;
                    --  
                );
            );
        );
    );
   
    --check:
    --print pairedRootRootList;
    apply(pairedRootRootList, rootlist->(assert (#rootlist==liftInfo#"requiredLatticeDimension"-1);));

    -- todo: liftAndLLLAndPairResult erklären!
    liftAndLLLAndPairResult := new MutableHashTable;
    for pos in 0..(#unknownList-1) do
        liftAndLLLAndPairResult#(unknownList#pos)=(minimalPolynomialsTable#(unknownList#pos)#0,minimalPolynomialsTable#(unknownList#pos)#1,pairedRootRootList#pos);
   
    roots := apply(#(liftAndLLLAndPairResult#(unknownList#0)#2),pos->apply(unknownList,unknown->liftAndLLLAndPairResult#unknown#2#pos));

    rootList := apply( roots,root->  matrix { root   } );
    
    return (minimalPolynomialsTable, rootList, liftInfo);
    -- todo: anonymer Rückgabetyp ist nicht gut. Allerdings gibt es erstmal ueberall Probleme, wenn Rückgabetyp hier geändert wird.
)


-- idealPointsApproxData: creates result data structure (a hash table)
--precondition: root coordinates corresponds to unknown list entries and unknownList is the same as gens ring equationsIdeal
idealPointsApproxData=( systemData, solutionPoint, minPolyTable, approxSolutions, mergedLiftInfo, unknownList )->
(
    approxSolutionData := new MutableHashTable;
    approxSolutionData#"systemData"  =  systemData ;
    approxSolutionData#"solutionPoint"     =  solutionPoint ;

    approxSolutionData#"approxVanishingSetElems"  =  approxSolutions ;
    approxSolutionData#"minPolyTable"     =  minPolyTable ;
    approxSolutionData#"mergedLiftInfo"  =  mergedLiftInfo ;
   
    
    errorList := {};
    bitPrecision := precision approxSolutions#0;

    complexRing := CC_bitPrecision[gens ring unknownList#0 ];

    for unknown in keys minPolyTable do
    (
        for root in approxSolutions do
        (
            residue := substitute(  minPolyTable#unknown#1,    root );
            errorList = append(errorList, abs( sub(residue, coefficientRing complexRing))   );
        );
    );
    approxSolutionData#"maxResidue" = max( errorList );
    
    approxSolutionData#"dataType" = "IdealPointsApprox";
    return new HashTable from approxSolutionData;
)

approxComplexSolutions = method (Options=>{"options"=>new LiftOptions,"decimalPrecision"=>10});

approxComplexSolutions (Ideal, Matrix) := opts->(inputIdeal, solutionPoint)->
(   
	return approxComplexSolutions( blackBoxIdeal(inputIdeal), solutionPoint );
)

-- may have problems in case unknownList is different from (gens ring inputIdeal order). Therefore removed the parameter for a moment.
approxComplexSolutions (HashTable, Matrix) := opts->( systemData, solutionPoint)->
(   

    --systemData.unknowns;

    --unknownList := gens ring inputIdeal;
    --unknowns := null;
    --unknowns = getGlobalSymbol "unknowns";
    --unknownList := systemData#unknowns;
    
    unknownList := systemData.unknowns;

    if (not pariGpIsPresent() ) then 
        error "please install Pari/GP or build Macaulay2 with pari (add 'pari' to the --enable-build-libraries configure parameter. ) to use 'approxComplexSolutions'";

    --localDecimalPrecision := opts#"options"#"decimalPrecision";
    localDecimalPrecision := opts#"decimalPrecision";

    rootCalculator := opts#"options"#"rootCalculator";

    characteristic:=char ring solutionPoint;


    (minimalPolynomialsTable, mergedLiftInfo ) := computeMinPolys( systemData, solutionPoint, unknownList, 
                                    "options" => opts#"options" );

    pairedRootRootList := new MutableList;

    opts#"options"#"initialLiftDepth"   = mergedLiftInfo#"maxLiftDepth"+1 ; --# is a heuristic. could be suboptimal for generic problems.
    opts#"options"#"initialLatticeDim"  = mergedLiftInfo#"requiredLatticeDimension" ;

    -- well, what should we do, return some info or return null in case something failed?
    if (mergedLiftInfo#"requiredLatticeDimension"===null) then  
    (
        return null;
        --return idealPointsApproxData( systemData, solutionPoint, minimalPolynomialsTable, {}  , mergedLiftInfo, unknownList );
    );
      
      --opts.logger(1, "------------------------pairing part---------------------------") ;
    
    rootListList := apply(unknownList, unknown-> rootCalculator(minimalPolynomialsTable#unknown#1, localDecimalPrecision) );   

    --for rootList in  rootListList do  (assert (#rootList==liftInfo#"requiredLatticeDimension"-1); );

    pairedRootRootList#0 = rootListList#0;
    operationRootList:={};

    operationAdd := (a,b)->(return a+b;);
    operationSub := (a,b)->(return a-b;);
    operationRnd := (a,b)->(return a+(random(ZZ/characteristic))*b;);
    --operationInputList:={operationAdd,operationRnd,operationSub};
    operationInputList := apply(characteristic-1, num-> ((a,b)->(return a+((num+1)_(ring b))*b;)));

    operationUsedList:={ };


    referenceRoots := rootListList#0;

    preApproxSolutions := new MutableList from apply(referenceRoots , n-> {{ n}} );

    unknown := unknownList#0;

    for unknownIdx in 1..#unknownList-1 do
    (
        currentCoordinatePaired := false;

         -- todo: check if pairing succeeded for each unknownIdx 
        for operation in operationInputList do
        (
            newUnknown := operation ( unknown , unknownList#unknownIdx );
            --opts.logger(2, Concatenation("newUnknown: ", String(newUnknown ) ) ) ;

            opts#"options"#"initialLatticeDim"=  1+ #preApproxSolutions ;
            
            -- adjust 'maxLatticeDim': the worst situtation would be if each root in preApproxSolutions
            --                          is compatible with each root in 'rootListList#unknownIdx'

            opts#"options"#"maxLatticeDim" =  1+ (#preApproxSolutions)*( #(rootListList#unknownIdx) ) ;
          
            --opts.logger(1, Concatenation("opts.maxLatticeDim: ", String(opts.maxLatticeDim ) ) ) ;
           (compatibilityMinimalPolynomials, compatLiftInfo ) := computeMinPolys( systemData, solutionPoint, {newUnknown}, 
                            "options" => opts#"options");
            
           if (compatLiftInfo#"requiredLatticeDimension"===null) then continue; --failed
        
            
            --opts.logger(1, Concatenation(" ----------------pairing variable ",String(unknownIdx) ) );
          
            unknownRoots := rootCalculator( compatibilityMinimalPolynomials#newUnknown#1, localDecimalPrecision);
                                             

            compMatrix := computeRootCompatibility( referenceRoots,rootListList#unknownIdx , unknownRoots  ,operation, opts#"options"#"maxPairingTolerance" );
            --print "compMatrix: "; print (toExternalString(compMatrix));
            if compMatrix===null then continue;  -- do not have a test for this!

           
            --opts.logger(2, "---------------------------compatibility matrix---------------------------------");
            --opts.logger(2, String(compMatrix) );
            tmppreApproxSolutions := new MutableList from apply(#unknownRoots, i->{}); --List([ 1..Size(unknownRoots)], n->[] );
                
            for row in 0..(numRows    compMatrix)-1 do
            for col in 0..(numColumns compMatrix)-1 do
            (
                if compMatrix_(row,col)>0 then
                    for entry in preApproxSolutions#row do
                    (
                        entryCopy :=  copy new List from entry; -- probably good enough here, but would be nice to have deepCopy()
                        --entryCopy :=  value toExternalString new List from entry;
                        entryCopy= append( entryCopy, rootListList#unknownIdx#col );
                        tmppreApproxSolutions#( compMatrix_(row,col) -1 ) = append( tmppreApproxSolutions#( compMatrix_(row,col)-1 ),  entryCopy  );
                    );
            );
            preApproxSolutions = tmppreApproxSolutions;
            --operationUsedList = append( operationUsedList, {operation} );
            referenceRoots = unknownRoots;
            unknown = newUnknown;
            currentCoordinatePaired = true;
            mergedLiftInfo = mergeLiftInfo( compatLiftInfo , mergedLiftInfo);
            --opts.logger(1, " ----------------pairing success---------------------------\n");
            --opts.logger(1, Concatenation("unknownIdx: ", String(unknownIdx) ) );
     
            break;        
           
        );-- end for
        if not currentCoordinatePaired then
        (
                error "current coordinate not paired";
                --opts.logger(0, Concatenation("pairing failed for indeterminate ", String(unknownIdx) ));
                return null;
        );
    );-- end for
    --opts.logger (1, " ---------------- All variables paired !---------------------------\n");
    --# todo: save input parameters in the result or not?
    --# debugInfo := Immutable (rec ( operationsUsedForPairing := operationUsedList ) );
   

    preApproxSolutions = flatten new List from preApproxSolutions;

    preApproxSolutions    = apply(preApproxSolutions, i-> matrix{i} );
    
    return idealPointsApproxData(  systemData, solutionPoint, minimalPolynomialsTable, preApproxSolutions  , mergedLiftInfo, unknownList );
    
    --return approxComplexSolutionData;
)


doc ///
    Key
        approxComplexSolutions
        (approxComplexSolutions ,Ideal, Matrix)
        (approxComplexSolutions ,HashTable, Matrix)
    Headline
        given a solution over a prime field for a system of equations, compute the corresponding complex solutions
    Inputs
         blackBoxIdeal: HashTable 
            an ideal blackbox (see @TO BlackBoxIdeals@ ) or the equations ideal ( only integer coefficient ring is supported ) 
        solutionModPrime: Matrix
            an element of the ideal vanishing set over a prime field
        "options"=> LiftOptions
            . If a computation takes too long, consider to customize {\tt options}.
        "decimalPrecision"=> ZZ
            . Decimal precision for the polynomial root computation step.
--        unknownList: List
--            the unknowns of interest (element of 'gens ring equationsIdeal')
    Outputs
        : HashTable
            contains complex solutions {\tt 'approxVanishingSetElems'} and intermediate result {\tt 'minPolyTable'}, 
            the minimal polynomials for the unknowns
    Usage
         complexSolutionsData = approxComplexSolutions( systemData, solutionModPrime, unknowns )
         complexSolutionsData = approxComplexSolutions( systemData, solutionModPrime, unknowns,"options"=>options, "decimalPrecision"=> decimalPrecision )

    Description
        Text
            Given a solution  over a prime field for an equation system, compute the corresponding complex solutions approximation.\break
            The function uses @TO computeMinPolys@ , @TO computeRootsWithGP@ and  @TO computeRootCompatibility@ . \break
            See also the article ... with the algorithm description on arXiv.    
        Example          
        Text
           \break  Example: find (subset of the) vanishing set for an ideal (here IFZ)
        Example          
            RQ = QQ[x,y];
            FQ = { x-1/3, y-1/5 };          
            IFQ = ideal FQ;        
            equationsIdeal =  clearCoeffDenominators IFQ
        Text
           \break the solutions over a finite field can be found via "brute force" -  omitted here
        Example          
            K = ZZ/11;
            solutionModPrime = matrix{ {4_K , -2_K} };  
            assert( sub(equationsIdeal, solutionModPrime) ==0 );   
        Text    
            \break compute the corresponding complex approximate solutions
        Example          
            complexSolutionsData = approxComplexSolutions( equationsIdeal, solutionModPrime, "decimalPrecision"=>20 );
            peek complexSolutionsData
            rootList =complexSolutionsData#"approxVanishingSetElems";
        Text
            \break check result correctness
        Example          
            complexRing := ring rootList#0[ gens ring IFQ ] ;
            sub( equationsIdeal, sub( rootList#0, complexRing ))
    Caveat
            as ideal coefficient ring currently only integers (ZZ) are supported and the solution set has to be 0-dimensional. 
///

-- todo: needs a test with more than one solution and one pairing - because I do not trust 'copy'
testApproxComplexSolutions = ()->
(
    x := null; x = symbol x;
    y := null; y = symbol y;
    prime := 7;

    RQ := QQ[ x,y ];
    x = (gens(RQ))#0;
    FQ := { (1/3*x+1/2), ((y-1/2) ) } ; 
    IFQ := ideal FQ;
    -- reformulate with integer coeff ring to get correct results:
    IFZ := clearCoeffDenominators(IFQ);  

    x = (gens ring IFZ )#0;
    y = (gens ring IFZ )#1;
    solutionQQ := matrix{{-3/2,1/2}} ;
    solution := sub( solutionQQ ,ZZ/prime) ;


    -- root gluing test:
    result := approxComplexSolutions( IFZ, solution  );
)

TEST ///
debug padicLift
padicLiftProtect()
testApproxComplexSolutions()
///
