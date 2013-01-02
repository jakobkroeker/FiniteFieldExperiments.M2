-- calculate p-adic lifts
--

--Questions: (ask in Forum): how to substitute a (nested) QQ coefficient ring to ZZ?
-- remark: for hashtable keys mainly strings are used, because string usage prevents introducing bugs 
--         caused by shadowing or overwriting symbols by local variables, etc.

-- todo: black box functionality!
-- replace assert messages with error messages.
-- todo: write lift tests, where the inteterminate order is different from ring variable order.
-- todo: inkonsequent: hashtables mit strings als Schlüssel und Gleichzeitig hashtables mit Symbolen als Schlüssel
-- todo: document tests
-- todo: clean design for message printing.
-- todo: report the 'factor 0'-bug, if still present.
-- todo: a test for failing lift situation 

-- checklist
--    indent: no Tabs, indent-width: 4.


newPackage(
     "padicLift",
     Version => "0.9.5", 
     Date => "05.11.2011",
     Authors => {{
           Name => "Jakob Kroeker", 
           Email => "kroeker@uni-math.gwdg.de", 
           HomePage => "http://www.crcg.de/wiki/User:Kroeker"},{
           Name => "Hans-Christian Graf v. Bothmer", 
           Email => "bothmer@math.uni-hannover.de", 
           HomePage => "http://www.crcg.de/wiki/Bothmer"}    
      },
     Configuration => {"gppath" =>  "gp" },
     Headline => "p-adic lift package",
     DebuggingMode => true,
     AuxiliaryFiles=>true
)


export{
    gppath,
    ReducedPadicLiftResult,
    LiftOptions,
    liftPoint,
    nextLift,
    computeMinPolys,
    approxComplexSolutions,
    approxComplexSolutionsOld,
    computeSingleMinPoly,
    plusOne,
    checkLiftOptions,
    computeRootsWithGP,
    createLiftOptions,
    pariGpIsPresent,
    disposeRationalCoeffs
}
    
 --protect the symbols for checking package correctness: no symbol variables should be overwritten by accident!
-- all in the package used hash table keys should be contained in this list.
-- e.g. in a  statement like 'liftOption.maxLiftDepth' maxLiftDepth should be a symbol variable without assigned value!
padicLiftProtect = ()->
(
    protect GlobalInternalPadicLiftResultVariable;
    protect unchanged;
    protect normalized;
    protect norms;
    protect normalizedNorms;
    protect initialLiftDepth;
    protect maxLiftDepth;
    protect initialLatticeDim;
    protect latticeDimIncrementFkt;
    protect maxLatticeDim;
    protect verbose;
    protect minColumnNormDistanceFactor;
    protect foundMinPolynomialCandidate;
    protect reducedLatticeBasis;
    protect latticeBasisVectorNormsList;
    protect currentLatticeDim;
    protect liftInfo;
    protect polynomial;
    protect unknown;
    protect maxLatticeDimension;
    protect requiredLatticeDimension;
    protect reductionOpts;
    protect liftOptions;
 
    protect tolerance;
    protect decimalPrecision;
)

padicLiftExport = ()->
(
    protect GlobalInternalPadicLiftResultVariable;
    exportMutable( unchanged);
    exportMutable( normalized);
    exportMutable( norms);
    exportMutable( normalizedNorms);
    exportMutable( initialLiftDepth);
    exportMutable( maxLiftDepth);
    exportMutable( initialLatticeDim);
    exportMutable( latticeDimIncrementFkt);
    exportMutable( maxLatticeDim);
    exportMutable( verbose);
    exportMutable( minColumnNormDistanceFactor);
    exportMutable( foundMinPolynomialCandidate);
    exportMutable( reducedLatticeBasis);
    exportMutable( latticeBasisVectorNormsList);
    exportMutable( currentLatticeDim);
    exportMutable( liftInfo);
    exportMutable( polynomial);
    exportMutable( unknown);
    exportMutable( maxLatticeDimension);
    exportMutable( requiredLatticeDimension);
    exportMutable( reductionOpts);
    exportMutable( liftOptions);
    exportMutable( tolerance);
    exportMutable( decimalPrecision);  

)

undocumented { 
    approxComplexSolutionsOld,
    unchanged,
    normalized,
    initialLiftDepth,
    maxLiftDepth,
    initialLatticeDim,
    latticeDimIncrementFkt,
    maxLatticeDim,
    verbose,
    minColumnNormDistanceFactor,
    foundMinPolynomialCandidate,
    reducedLatticeBasis,
    latticeBasisVectorNormsList,
    currentLatticeDim,
    liftInfo,
    polynomial,
    unknown,
    maxLatticeDimension,
    latticeBasisVectorToPolynomial,
    configurePari,
    decimalPrecision,
    liftOptions,
    normalizedNorms,
    norms,
    reductionOpts ,
    requiredLatticeDimension,
    tolerance,
    createLiftOptions
}

-- swith between protect and export - both are not possible!

--padicLiftProtect() -- protect the symbols for checking package correctness: no symbol variables should be overwritten by accident!
padicLiftExport(); -- export the symbols to make the package work 

-- a package global  for indeterminate usage ...hopefully not too bad design...?

GlobalInternalPadicLiftResultVariable := null;
GlobalInternalPadicLiftResultVariable = symbol GlobalInternalPadicLiftResultVariable;
GlobalInternalPadicLiftResultVariable = (gens (ZZ[GlobalInternalPadicLiftResultVariable]))#0;



beginDocumentation()


load "./padicLift/univarPolRoots.m2";
load "./padicLift/rootPairing.m2";



doc ///
    Key
        padicLift
    Headline
          lifting of polynomial system solutions over a prime field to an extension ring of rational numbers.
    Description
        Text
            Implements lifting of a smooth isolated polynomial system solution over a finite field to an extension of rationals \break \break
            
             \bullet \, for simple Hensel lifting see @TO liftPoint@ and  @TO nextLift@  \break
             \bullet \, for minimal polynomial computation given a finite field solution for an equation system see  @TO computeMinPolys@ \break
             \bullet \, the roots of the minimal polynomials above can be approximated by @TO approxComplexSolutions@ \break
             \bullet \, the package also provides  @TO computeRootsWithGP@, 
                        an interface to Pari for computing (complex) roots of univariate polynomials  .
    Caveat
        Does not implement yet black box interface.
///





-- removes a constant factor from a product - returns the new value. 
-- Attention: factoring the new value may result in  a Product with constant factors again depending on the factoring algorithm.
removeConstantFactors  = method();

removeConstantFactors(RingElement) := Product=> (ringElement)->
(
    result := null;
    localFactors := factor ringElement;
    return value removeConstantFactors(localFactors);
)

-- removes a constant factor from a Product.  - returns left factors as a Product
removeConstantFactors(Product) := Product=> (factors)->
(
    result := null;
    localFactors := factor value factors;
    for currentfactor in localFactors do
    (
        if ((degree currentfactor#0)#0!=0 ) then 
        (
            if (result===null) then 
            result=Product(currentfactor)
            else
            result=result*currentfactor;
        );
    );
    return result;
)


testRemoveConstantFactors=()->
(
    x := null; x=symbol x;
    rng := ZZ[x];
    x = (gens rng )#0;
    polynomial := 3*(x-1)*(x+2);
    monicPol := removeConstantFactors(factor polynomial);
    assert( value monicPol==(x-1)*(x+2) );
    assert(#monicPol==2);
    --------
    monicPol = removeConstantFactors(factor polynomial);
    assert( value monicPol==(x-1)*(x+2) );
    assert(#monicPol==2);
    --------
    polynomial  = 3 *(x+2);
    monicPol  = removeConstantFactors(factor polynomial);
    assert( value monicPol== (x+2) );
    assert(#monicPol==1);
    --------
    polynomial  = Product(3 *(x+2));
    monicPol  = removeConstantFactors(polynomial);
    assert( value monicPol== (x+2) );
    assert(#monicPol==1);
)

-- 'liftStep':   lifts an element of the vanishing set of the ideal sub(equationsIdeal,  ZZ[]/currchar)
-- to ZZ[]/currchar*currchar with currchar=char ring solution; see Hensel lifting for p-adic numbers. 
-- The default case is that the jacobian(ideal) is square and has full rank. 
-- Other cases are not considered in full extent, even if Hensel lifting would work for a specific case in theory.
-- The computation expects an ideal (polynomials) with integer coefficients; for ideals with rational coefficients 
-- use disposeRationalCoeffs() first!
-- Parameters:
-- equationsIdeal: equations ideal with integer coefficient ring ,
-- TransposedJacobian ( in Macaulay the jacobian is always transposed)
-- 'vanishingCoordinates' : element of a vanishing set of ideal mod prime^k  -  'sub( equationsIdeal, vanishingCoordinates )' should be zero.
--Returns:
-- 'nextVanishingCoordinates': element of a vanishing set of ideal mod (prime^k)^2  =>  'sub( equationsIdeal, nextVanishingCoordinates )' will be zero.
-- 
liftStep = (equationsIdeal, TransposedJacobian, vanishingCoordinates ) ->
(
    assert( sub( equationsIdeal, vanishingCoordinates)==0 );

    if  ((coefficientRing ring TransposedJacobian ) =!= ZZ) then error " liftStep works correctly only for ideals in ZZ";
    currchar := char ring vanishingCoordinates;
    nextchar := currchar*currchar;
    nextLiftDestRing := ZZ[]/nextchar; 
    localVanishingCoordinates := sub( vanishingCoordinates, nextLiftDestRing );
    -- TransposedJacobianAtSolution := sub( TransposedJacobian , localVanishingCoordinates);
    TransposedJacobianAtSolution  := sub( TransposedJacobian , vanishingCoordinates); -- this is sufficient.

    prime := ( factor currchar)#0#0;
     
    rightHandSide := transpose gens sub( equationsIdeal, localVanishingCoordinates );
    JacobianInverseAtSolution := null;
        
    if (rank source TransposedJacobianAtSolution ) > (rank image sub(TransposedJacobianAtSolution,ZZ) ) then
    (
        -- compute the left inverse
        JacobianInverseAtSolution =  transpose (TransposedJacobianAtSolution)^-1;        
        assert(JacobianInverseAtSolution*(transpose TransposedJacobianAtSolution)==1);
    );

    if (rank source TransposedJacobianAtSolution) <= (rank image sub(TransposedJacobianAtSolution,ZZ) ) then
    (
        -- compute the  right inverse
        JacobianInverseAtSolution = (transpose TransposedJacobianAtSolution)^-1;
        assert( (transpose TransposedJacobianAtSolution)*JacobianInverseAtSolution == 1 );
    );
           
    --firstCorrectionPart  := -transpose (JacobianInverseAtSolution* rightHandSide );
    firstCorrectionPart  := -transpose (sub(JacobianInverseAtSolution,nextLiftDestRing)* rightHandSide );

    --JacobianKernelAtSolution  := syz transpose TransposedJacobianAtSolution;
    --kerdim := rank source JacobianKernelAtSolution;
    --randomvector := matrix { apply(kerdim, i->(random(char nextLiftDestField))_nextLiftDestField ) };
    --randomvector =sub (randomvector, nextLiftDestField);
    --secondCorrectionPart := currchar_nextLiftDestField*randomvector*(transpose JacobianKernelAtSolution);
    
    nextVanishingCoordinates :=  firstCorrectionPart  + localVanishingCoordinates;
    assert( sub( equationsIdeal, nextVanishingCoordinates)==0 );

    return   nextVanishingCoordinates;
)


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


-- disposeRationalCoeffs converts an ideal with rational coefficients to an ideal with integer coefficients while preserving the vanishing set.
-- e.g. if sub(IdealWithRationalCoeffs,point)==0, then  sub( disposeRationalCoeffs(IdealWithRationalCoeffs),point)==0 and vice versa
disposeRationalCoeffs = method();

disposeRationalCoeffs (Ideal)  :=  Ideal =>  (IdealWithRationalCoeffs)->
(
    if (coefficientRing ring IdealWithRationalCoeffs=!=ZZ and coefficientRing ring IdealWithRationalCoeffs=!=QQ) then
    error("expected rationals as coefficient ring!");
    dstrng := ZZ[gens ring IdealWithRationalCoeffs];
    modgens := apply(flatten entries gens IdealWithRationalCoeffs, i->polynomialLCMDenominator(i)*i );
    return sub(ideal modgens,dstrng );
)

doc ///
    Key
        disposeRationalCoeffs        
        (disposeRationalCoeffs, Ideal )
    
    Headline
        convert an ideal with rational coefficients to an ideal with integer coefficients
    Usage   
        disposeRationalCoeffs(IdealInQQ)
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
            IFZ = disposeRationalCoeffs(IFQ)
///

testDisposeRationalCoeffs =()->
(
    x := null;  x=symbol x;
    y := null;  y=symbol y;
    RQ := QQ[x,y];
    x = (gens(RQ))#0;
    FQ := { (1/3*x+1) ,  (y+1/2)}; 
    IFQ := ideal FQ;
    IFZ := disposeRationalCoeffs(IFQ);  
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

testTensoredDisposeRationalCoeffs =()->
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
    IFZ := disposeRationalCoeffs(IFQ);  
    FZ := (entries (gens IFZ)_0)#0;
    assert(   FZ == sub(x*z+3*z,ring FZ)   ) ;       
)


-- 
testLiftStep = ()->
(
    x:=null; x=symbol x;
    y:=null;  y=symbol y;
    RZ := ZZ[x,y];
    FZ := 33*x^3+19*x^2-81*x-4;          
    IFZ := ideal FZ;
        
    K := ZZ/11;
    solution1 := matrix{ {1_K , 1_K} }; 
    solution2 := matrix{ {5_K, 2_K} };

    assert( sub(IFZ, solution1) ==0 );   
    assert( sub(IFZ, solution2 ) ==0);        

    nextApprox := nextLift(IFZ,solution2);
    assert (sub (IFZ, nextApprox)==0);       

    nonsolution := matrix{ {2_K , 2_K} }; 
    assert(sub(IFZ, nonsolution)!=0);

-- test different variable order

    RZ = ZZ[y,x];
    FZ = 33*x^3+19*x^2-81*x-4;          
    IFZ = ideal FZ;
        
    K = ZZ/11;
    solution1 = matrix{ {1_K , 1_K} }; 
    solution2 = matrix{ {5_K, 2_K} };
    solution2 = matrix{ {2_K, 5_K} };

    assert( sub(IFZ, solution1) ==0 );   
    assert( sub(IFZ, solution2 ) ==0);        

    nextApprox = nextLift(IFZ,solution2);
    assert (sub (IFZ, nextApprox)==0);       

)

-- 'nextLift()': interface;  lifts element of ideal IZ in ZZ/prevchar to ZZ/nextchar
--preconditions : ideal IZ is in ZZ  ans solution is element of the correct ring (  'sub(IZ,solution )==0' !)
--
-- Parameters:
-- IZ=ideal in ZZ, 
-- JZ: jacobian(IZ) ,
-- solution: current lift as element of ZZ/(p^k)
--
nextLift =method();

nextLift (Ideal,Matrix, Matrix)  := Matrix=> (equationsIdeal, TransposedJacobianOfEquations , vanishingCoordinates) -> (
    assert( sub( equationsIdeal,vanishingCoordinates ) == 0 );
    assert( (char ring vanishingCoordinates ) != 0 );
    nextVanishingCoordinates := liftStep(equationsIdeal, TransposedJacobianOfEquations, vanishingCoordinates);
    --print(nextVanishingCoordinates);
    assert( sub( equationsIdeal, nextVanishingCoordinates)==0 );
    return nextVanishingCoordinates;
)


-- todo: change second parameter type to vector?
nextLift (Ideal,Matrix) :=  Matrix=> (equationsIdeal,  vanishingCoordinates) -> (
    return nextLift(equationsIdeal,jacobian equationsIdeal, vanishingCoordinates);
)

doc ///
    Key
        nextLift        
        (nextLift, Ideal,Matrix)
        (nextLift, Ideal,Matrix,Matrix)
    Headline
        lift a solution point mod prime^k to  mod (prime^k)^2  
    Usage   
        nextLift(equationsIdeal,solutionPoint)
    Inputs
        equationsIdeal: Ideal
             ideal with integer coefficient ring.
        vanishingCoordinates: Matrix
             element of the ideal vanishing set over a finite field in matrix form
    Outputs
        : Matrix
            lifted solutionPoint in matrix form. 
            If the input { \tt vanishingCoordinates } are vanishing coordinates mod p, the result are vanishing coordinates mod p^2.
    Description
        Example          
        Text
           \break  Example: lift a finite field vanishingCoordinates for an ideal in ZZ (here IFZ) \break
           1. define the problem
        Example          
            RQ = QQ[x,y];
            FQ = 33/4*x^3+19/4*x^2-81/4*x-1;          
            IFQ = ideal FQ
        Text
            1a. to do lifting correctly in Macaulay2, we need to translate the problem to a ring with integer coefficients (ZZ):
        Example
            IFZ = disposeRationalCoeffs(IFQ)
        Text
           \break 2. the vanishingCoordinatess over a finite field can be found via brute force -  omitted here
        Example          
            K = ZZ/11;
            vanishingElem1 = matrix{ {1_K , 1_K} }; vanishingElem2 = matrix{ {5_K, 1_K} };
        Text
           \break 3. check solutions
        Example
            assert( sub(IFZ, vanishingElem1) ==0 );   assert( sub(IFZ, vanishingElem2 ) ==0);        
        Text    
            \break 4. compute next padic approximation
        Example
            nextApprox = nextLift(IFZ,vanishingElem2)
            sub (IFZ, nextApprox)
            assert (sub (IFZ, nextApprox)==0);                 
///

-- 'liftPoint()': interface;  lift element 'solution' of the vanishing set of the ideal IZ in ZZ/primepower to ZZ/(primepower^(2^numLiftDepth))
-- returns  vanishingCoordinates mod ZZ/(primepower^(2^numLiftDepth))
--
-- Parameters:
-- equationsIdeal  IdealWithIntegerCoeffRing
-- vanishingCoordinates:  element of the vanishing set of ideal IZ  mod a 'prime'; 
-- numLiftDepth: solution is lifted to  ZZ/(prime^(2^numLiftDepth))
--  
liftPoint = (equationsIdeal,  vanishingCoordinates , numLiftDepth) -> (
    --assert(isPrime char ring vanishingCoordinates);
    assert((char ring vanishingCoordinates) !=0 ); 
    -- if characteristic was zero, the next test (factor char...) could also crash Macaulay2 (factoring zero. bug?)
    assert(#(factor char ring vanishingCoordinates) == 1 ); -- characteristic of the solution is a power of the prime
    assert( sub (equationsIdeal, vanishingCoordinates ) == 0 );
    IdealTransposedJacobian := jacobian equationsIdeal;
    currLiftDepth := 0;
    localVanishingCoordinates := vanishingCoordinates;
    while currLiftDepth < numLiftDepth do
    (
        localVanishingCoordinates = nextLift( equationsIdeal, IdealTransposedJacobian, localVanishingCoordinates);
        --print(localVanishingCoordinates);
        currLiftDepth = currLiftDepth + 1;
    );
    return localVanishingCoordinates;
)
doc ///
    Key
        liftPoint
    Headline
          lift a vanishingElem of the  equationsIdeal ( mod p) to  mod p^{2^{numLiftDepth}}.
    Usage   
        liftPoint( equationsIdeal, vanishingElem, numLiftDepth)
    Inputs
        equationsIdeal: Ideal
        vanishingElem: Matrix
             element of the ideal vanishing set over a finite field in matrix form
        numLiftDepth: ZZ 
            If the input vanishingElem is an element of the ideal vanishing set mod p, 
            the result is a an element of the ideal vanishing set mod p^{2^{numLiftDepth}}.
    Outputs
        : Matrix
            lifted vanishingElem. If the input vanishingElem is a element of the ideal vanishing set mod p, 
            the result is a vanishingElem mod p^{2^{numLiftDepth}}.
    Description
        Example          
        Text
           \break  Example: lift a finite field vanishingElem of an ideal in QQ (here IFQ)
           1. formulate the problem 
        Example          
            RQ = QQ[x,y];
            FQ = 33/4*x^3+19/4*x^2-81/4*x-1;          
            IFQ = ideal FQ;
        Text
           \break 2. reformulate the equations in  IFQ without denominators - mandatory for {\tt liftPoint }
        Example          
            IFZ = disposeRationalCoeffs(IFQ)
            
        Text
           \break 3. the vanishingElems over a finite field can be found via "brute force" -  omitted here
        Example          
            K = ZZ/11;
            vanishingElem1 = matrix{ {1_K , 1_K} }; vanishingElem2 = matrix{ {5_K, 1_K} };
        Text
           \break 4. check that the points are elements of the ideal vanishing set mod 11:
        Example
            assert( sub(IFZ, vanishingElem1) ==0 );   assert( sub(IFZ, vanishingElem2 ) ==0);
        Text    
            \break 5. compute the padic approximation  mod prime^(2^4). 
        Example
                liftDepth := 4;
                liftResult = liftPoint( IFZ, vanishingElem2, liftDepth )
                sub (IFZ, liftResult)
                assert ( sub (IFZ, liftResult)==0 );                 
///


testLiftPoint = ()->
(
  --test nonstandard situation where the jacobian is not square.
    x:=null; x=symbol x;
    y:=null;  y=symbol y;
    RQ := QQ[x,y];
    x  = (gens RQ)#0;
    y  = (gens RQ)#1;
    FQ := 33/4*x^3+19/4*x^2-81/4*x-1;          
    IFQ := ideal FQ;
        
    IFZ := disposeRationalCoeffs(IFQ);
    prime := 11;
    K := ZZ/prime;
    --some solutions of the equations given by the generators of the ideak IFZ:
    solution1 := matrix{ {1_K , 1_K} }; 
    solution2 := matrix{ {5_K, 2_K} };

    assert( sub(IFZ, solution1) ==0 );   
    assert( sub(IFZ, solution2 ) ==0);        

    liftDepth:=10;
    nextApprox := liftPoint(IFZ,solution2 ,liftDepth);
    assert (prime^(2^liftDepth)==char ring nextApprox);

    assert (sub (IFZ, nextApprox)==0);       

    nonsolution := matrix{ {2_K , 2_K} }; 
    assert(sub(IFZ, nonsolution)!=0);

-- test different variable order

    RQ = QQ[y,x];
    x  = (gens RQ)#1;
    y  = (gens RQ)#0;
    FQ = 33/4*x^3+19/4*x^2-81/4*x-1;          
    IFQ = ideal FQ;
    IFZ = disposeRationalCoeffs(IFQ);  
    K = ZZ/prime;
    solution1 = matrix{ {1_K , 1_K} }; 
    solution2 = matrix{ {5_K, 2_K} };
    solution2 = matrix{ {2_K, 5_K} };

    assert( sub(IFZ, solution1) ==0 );   
    assert( sub(IFZ, solution2 ) ==0);        

    nextApprox = liftPoint(IFZ,solution2,liftDepth);
    assert (prime^(2^liftDepth)==char ring nextApprox);
    assert (sub (IFZ, nextApprox)==0);       

)


-- returns a list of the euclidean matrix column norms
computeColumnNorms=(M)->
(
    lM  := sub(M,ZZ);
    MM := transpose lM * lM;
    rs := rank source MM;
    --normlist := apply(rs,i-> MM_i_i );
    --normlist := apply(rs,i->round(power(MM_i_i,1/2)));
    normlist := apply(rs,i->power(MM_i_i,1/2) );
    return normlist;
)


-- returns for a matrix M a data table with 
-- .norms           : the matrix M column norms
-- .normalizedNorms :  matrix M column norms divided by the smallest occured norm.
-- .min             : smallest matrix column norm
-- .max             : biggest matrix column norm
computeNormalizedColumnNormsData=(M)->
(   
    columnNormlist := computeColumnNorms(M);
    minNorm := min(columnNormlist);
    normalizedColumnNormlist :=  apply(columnNormlist, columnnorm-> round sub(columnnorm/minNorm,RR));
    result := new MutableHashTable;
    result.norms = columnNormlist;
    result.normalizedNorms = normalizedColumnNormlist;
    result.min = minNorm;
    result.max = max(columnNormlist);
    return result;    
)

plusOne=(k)->
(
    assert( class(k)===ZZ or class(k)===InfiniteNumber);
    return k=k+1;
)


plusTen=(k)->
(
    assert( class(k)===ZZ or class(k)===InfiniteNumber);
    return k=k+10;
)

incrementByDoubling=(k)->
(
    assert( class(k)===ZZ or class(k)===InfiniteNumber);
    if (k<40) then
        return k=2*k
    else
    k = k+10;
    return k;
)




--- constructLLLInputFromLift: construct such a lattice basis (matrix columns), that the shortest element of its LLL-reduced basis
--- may contain coefficients of the minimal polynomial for the unknown.
---
--- Remark1: the  degree for the minimal polynomial is guessed as 'currentLatticeDim-1'
--- Remark2: liftResult_i is the current padic approximation for the unknown  if unknown is the i-th ring variable 
--           (get ring variables with 'gens ring unknown')
---
-- Example: for currentLatticeDim=3 the   columns will be {unknown,-1,0}, {unknown^2,0,-1}, { (char ring liftResult),0,0 }
-- the (small) reduced vector will be a linear combination of the constructed basis and therefore of the form 
-- (  a*unknown + b*unknown^2 + c*(char ring liftResult),   -a,   -b ) 

constructLLLInputFromLift = (unknown, liftResult, currentLatticeDim )->
(
    liftResultOverZZ := sub( liftResult, ZZ ); 
    M :=   matrix{ apply(currentLatticeDim, i->unknown^i) } ;
    M  =  sub( M ,  liftResultOverZZ ); -- substitute the current padic approximation for the unknown
    M = M |  matrix{{ char ring liftResult }};
    sM := (syz M)^{0..(currentLatticeDim-1)}; -- delete last row
    return  sub(sM,ZZ);
)

-- construct a lattice basis as above (constructLLLInputFromLift) without using the syz Macaulay2 command. 
constructLLLInputFromLiftWithoutSyz = (unknown, liftResult, currentLatticeDim )->
(
    liftResultOverZZ := sub( liftResult, ZZ ); 
    M :=  sub( matrix{ apply(currentLatticeDim, i->unknown^i) },  liftResultOverZZ );
    M = M |  matrix{{ char ring liftResult }};
    sM := mutableMatrix apply(currentLatticeDim+1, i->apply(currentLatticeDim,i->0));
    for idx in 0..currentLatticeDim-1 do
    (
        sM_(0,idx) = -M_(0,idx+1);
        sM_(idx+1,idx) = M_(0,0);
    );
    sM = matrix sM;
    result :=  (sM)^{0..(currentLatticeDim-1)};
    return result;
)

------------------------------------------------------------------------------------------

LiftOptions = new Type of MutableHashTable;


createValidLiftOptionsSpeciman  = () ->
(
    preLiftOptions := new MutableHashTable;  
    preLiftOptions#"initialLiftDepth" = null;
    preLiftOptions#"maxLiftDepth"       =  null;
    preLiftOptions#"initialLatticeDim"  =  null;
    preLiftOptions#"latticeDimIncrementFkt" = null;
    preLiftOptions#"maxLatticeDim" = null;
    preLiftOptions#"verbose"        =  null;
    preLiftOptions#"minColumnNormDistanceFactor"= null;
    preLiftOptions#"decimalPrecision" = null;
    preLiftOptions#"maxPairingTolerance" = null;
    preLiftOptions#"LLLStrategy" = null;
    preLiftOptions#"rootCalculator" = null;
    return  preLiftOptions;
)


-- check the consistency of the lift options - a common error is to set not intended keys , 
-- e.g. by using 'liftOpts#initialLatticeDim' instead of 'liftOpts#"initialLatticeDim'"
checkLiftOptions = method();
checkLiftOptions (Thing) := (liftOpts)->
(
    --print opts;
    --print liftOpts;        
    assert( class(liftOpts#"initialLiftDepth")===ZZ );
    assert( class(liftOpts#"initialLatticeDim")===ZZ );
    assert( class(liftOpts#"decimalPrecision")===ZZ );
    assert( class(liftOpts#"maxPairingTolerance")===RR  );
    assert ( liftOpts#"maxPairingTolerance" > 0 );
    assert ( liftOpts#"decimalPrecision" >= 1 );
    assert( class(liftOpts#"verbose")===Boolean );
    assert( class(liftOpts#"minColumnNormDistanceFactor")===ZZ );
    --
    assert( class(liftOpts#"maxLiftDepth")===ZZ or class(liftOpts#"maxLiftDepth")===InfiniteNumber   );
    assert( class(liftOpts#"maxLatticeDim")===ZZ or class(liftOpts#"maxLatticeDim")===InfiniteNumber   );
    --
    assert ( liftOpts#"initialLatticeDim" >= 1 );
    assert(  liftOpts#"maxLatticeDim" >= liftOpts#"initialLatticeDim");
    --
    validLiftOptions := createValidLiftOptionsSpeciman();
    for key in keys liftOpts do
    (
        if (not validLiftOptions#?key) then
            error ("LiftOptions contains invalid key " | toString key |". Did you incorrectly access modify entry? - All table keys are expected as strings! ");
    )
)

 

new LiftOptions from MutableHashTable :=(A,mht)->
(
    checkLiftOptions(mht);
    return mht;
)



-- createLiftOptions declaration; internal
createLiftOptions = method(Options=>{
                    "initialLiftDepth"=>0,
                    "maxLiftDepth"=>infinity,
                    "initialLatticeDim"=>1,
                    "latticeDimIncrementFkt"=>plusOne,
                    "maxLatticeDim"=>infinity,
                    "verbose"=>false,
                    "minColumnNormDistanceFactor"=>10,
                    "decimalPrecision"=>16,
                    "LLLStrategy"=>NTL,
                    "maxPairingTolerance"=>0.01,
                    "rootCalculator"=>computeRootsWithGP
                                }
);

-- createLiftOptions definition; internal
-- note: the dummy parameter was used here, since it is syntactically(?) more difficult (or not possible?)
-- to define a method without a parameter

createLiftOptions(Thing) := opts->(dummy)->
(
    preLiftOptions := new MutableHashTable;  
    preLiftOptions#"initialLiftDepth"  =  opts#"initialLiftDepth"; 
    preLiftOptions#"maxLiftDepth"       =  opts#"maxLiftDepth"; 
    preLiftOptions#"initialLatticeDim"  =  opts#"initialLatticeDim"; 
    preLiftOptions#"latticeDimIncrementFkt" =  opts#"latticeDimIncrementFkt"; 
    preLiftOptions#"maxLatticeDim"  =  opts#"maxLatticeDim"; 
    preLiftOptions#"verbose"        =  opts#"verbose"; 
    preLiftOptions#"minColumnNormDistanceFactor"= opts#"minColumnNormDistanceFactor"; 
    preLiftOptions#"decimalPrecision" = opts#"decimalPrecision";
    preLiftOptions#"LLLStrategy"      = opts#"LLLStrategy";  
    preLiftOptions#"maxPairingTolerance"      = opts#"maxPairingTolerance";  
    preLiftOptions#"rootCalculator"   = opts#"rootCalculator";  
    return new LiftOptions from preLiftOptions;
)

new LiftOptions := (ObjType)->
(
    return createLiftOptions( null );
)


doc ///
    Key
        checkLiftOptions, Thing
    Headline
        check options consistency 
    Description
        Text
            checks consistency of @TO LiftOptions@
               
    Caveat
///


doc ///
    Key
        LiftOptions
    Headline
        options for p-adic lifting and LLL reduction 
    Usage
        liftOptions = new LiftOptions
    Description

        Text 
                Used in :   @TO computeMinPolys@,   @TO approxComplexSolutions@\break
                LiftOption keys are:  \break
                \,\, \bullet \, { \tt "initialLiftDepth" } (default: 0), an @TO2 {ZZ,"integer"}@ ,  a point from an ideal vanishing set mod prime is lifted to prime^{2^{initialLiftDepth}}, bevore a first lattice reduction attempt is started   \break
                \,\, \bullet \, { \tt "maxLiftDepth"} (default: infinity), an @TO2 {ZZ,"integer"}@ ,   if the minimal polynomial computation fails for a point lifted to prime^{2^{maxLiftDepth}}, the computation stops    @BR{}@
                \,\, \bullet \, { \tt "initialLatticeDim"} (default: 1), an @TO2 {ZZ,"integer"}@ ,      smallest value for { \tt latticeDimension } used in the lattice basis reduction step     \break
                \,\, \bullet \, { \tt "maxLatticeDim"} (default: infinity),  an @TO2 {ZZ,"integer"}@ ,      max value for latticeDim used in the lattice basis reduction step. \break
                                                                                       \,\, \,\, \,\, \,\, lattice reduction is performed until  minimal polynomial computation succeeded or until {\tt maxLatticeDim} is reached for the current (fixed) padic point approximation    \break
                \,\, \bullet \, { \tt "latticeDimIncrementFkt"} (default: plusOne), a @TO2 {MethodFunction,"method"}@ ,      Increment function for the {\tt latticeDimension } in the  lattice reduction step. Takes and returns an  @TO2 {ZZ,"integer"}@    \break
                \,\, \bullet \, { \tt "minColumnNormDistanceFactor"} (default: 10) , an @TO2 {ZZ,"integer"}@ ,      additional solution filter heuristic.  \break Expects that the the reduced lattice basis vector norms differ at least by a factor  {\tt minColumnNormDistanceFactor}.  Used to eliminate false positives.    \break
                \,\, \bullet \, { \tt "LLLStrategy" } (default: NTL), see @TO  {LLL} @ ,    the lattice reduction algorithm version     \break      
            \,\, \bullet \, { \tt "rootCalculator" } (default: @TO computeRootsWithGP@), an univariate polynomial root calculator      \break     
                \,\, \bullet \, { \tt "maxPairingTolerance" } (default: 0.01), maximal allowed error for root coordinate pairing.  (used in @TO approxComplexSolutions@ )    \break     
                \,\, \bullet \, { \tt "verbose" } (default: false), a @TO2 {Boolean,"boolean value"}@ ,  if {\tt true}, some  intermediate information is printed    \break   

        Example
        Text
            Usage example:
            
        Example
             liftOptions := new LiftOptions;
             peek  new LiftOptions
        Text
            \break   Set an option (here "maxLatticeDim") manually:
        Example
             liftOptions#"maxLatticeDim"= 79;
             checkLiftOptions liftOptions --ok
            --shape := computeShape(pol)
            --factor pol
    Caveat
        no setter and getter functions implemented. 
///


LatticeReductionResult = new Type of MutableHashTable;


new LatticeReductionResult from  HashTable := (ObjType,mht)->
(
    assert(mht.?foundMinPolynomialCandidate);
    assert(instance(mht.foundMinPolynomialCandidate,Boolean));
    
    assert(mht.?reducedLatticeBasis);
    assert(mht.reducedLatticeBasis===null or instance(mht.reducedLatticeBasis,Matrix));

    assert(mht.?latticeBasisVectorNormsList);
    assert(instance(mht.latticeBasisVectorNormsList,List));

    assert(mht.?currentLatticeDim);
    assert(instance(mht.currentLatticeDim,ZZ));
 
    return mht;    
)

doc ///
    Key
        LatticeReductionResult
    Headline
        lattice reduction result type 
    Description

        Text 
                LatticeReductionResult keys are:  \break \break
                \,\, \bullet \, { \tt "foundMinPolynomialCandidate" }, an @TO2 {Boolean,"boolean value"}@ ,  signals if a result candidate was found or not   \break
                \,\, \bullet \, { \tt "reducedLatticeBasis"}, an @TO Matrix @ ,  resulting @TO LLL@ matrix   \break
                \,\, \bullet \, { \tt "latticeBasisVectorNormsList"}, an @TO Matrix @ ,     list of the lll matrix column  norms     \break
                \,\, \bullet \, { \tt "currentLatticeDim"}, an @TO2 {ZZ,"integer"}@ ,   the used lattice dimension for reduction \break
    Caveat
///

---------------------------------------------------------

-- todo: wozu gehoert latticeBasisVectorNormsList, zum Ergebnis oder zu LiftInfo?

ReducedPadicLiftResult = new Type of MutableHashTable;

new ReducedPadicLiftResult from  HashTable := (ObjType,mht)->
(
    assert(mht.?foundMinPolynomialCandidate);
    assert(instance(mht.foundMinPolynomialCandidate,Boolean));
    
    --assert(mht.?reducedLatticeBasis);
    --assert(instance(mht.reducedLatticeBasis,Matrix));

    --assert(mht.?latticeBasisVectorNormsList);
    --assert(instance(mht.latticeBasisVectorNormsList,List));

    assert(mht.?liftInfo);
    assert(instance(mht.liftInfo,LiftInfo));

    assert(mht.?polynomial);
    assert(instance(mht.polynomial,RingElement));

    assert(mht.?unknown);

    return mht;    
)
doc ///
    Key
        ReducedPadicLiftResult
    Headline
        single variable lift result result type 
    Description
        Text 
                single variable lift result result type, see @TO computeSingleMinPoly@  \break
                ReducedPadicLiftResult keys are (please access values via '.' operator - accessing via '#' is unsafe !):  \break \break
                \,\, \bullet \, { \tt foundMinPolynomialCandidate }, an @TO2 {Boolean,"boolean value"}@ ,  signals if a result candidate was found or not   \break
                \,\, \bullet \, { \tt reducedLatticeBasis}, an @TO Matrix @ ,  resulting @TO LLL@ matrix   \break
                \,\, \bullet \, { \tt latticeBasisVectorNormsList}, an @TO List @ ,     list of the {\tt reducedLatticeBasis } column  norms     \break
                \,\, \bullet \, { \tt liftInfo}, an @TO LiftInfo@ ,   the required minimal lattice dimension and minimal lift depth \break
                \,\, \bullet \, { \tt polynomial}, an @TO RingElement@ ,  the roots of the polynomial are the solutions for the lifted unknown \break
                \,\, \bullet \, { \tt unknown}, an @TO RingElement@ ,  the lifted unknown variable \break
    Caveat
///

-- todo: intruduce liftResultSet and polynomials? - no , because LLLOutput may cost a lot of memory.

LiftInfo = new Type of HashTable;

new LiftInfo from HashTable := (A,ht)->
(
    assert(ht#?"maxLiftDepth");
    assert(ht#?"maxLatticeDimension");
    assert(ht#?"requiredLatticeDimension");

    assert(class ht#"maxLiftDepth"===ZZ);
    assert(class ht#"maxLatticeDimension"===ZZ);
    assert(class ht#"requiredLatticeDimension"===ZZ or   ht#"requiredLatticeDimension"===null);
  
    return ht;    
)

doc ///
    Key
        LiftInfo
    Headline
        lift info
    Description

        Text 
                LiftInfo keys are (please access via '.' operator - accessing via '#' is potentially unsafe!):  \break \break
                \,\, \bullet \, { \tt maxLiftDepth }, an @TO2 {ZZ,"integer"}@ ,  maximal used lift depth during computation   \break
                \,\, \bullet \, { \tt maxLatticeDimension}, an  @TO2 {ZZ,"integer"}@,  maximal used lattice dimension during lattice reduction \break
                \,\, \bullet \, { \tt requiredLatticeDimension}, an  @TO2 {ZZ,"integer"}@,  minimal required lattice dimension for successfull lattice reduction; null if unknown.  \break                                
    Caveat
///


createLiftInfo = method();
createLiftInfo := (maxLiftDepth,maxLatticeDimension,requiredLatticeDimension)->
(
    ht := new MutableHashTable;
    ht#"maxLiftDepth" = maxLiftDepth;
    ht#"maxLatticeDimension" = maxLatticeDimension;
    ht#"requiredLatticeDimension" = requiredLatticeDimension;
    return new LiftInfo from ht;
)



-- merging two lift statistics (usually for two different unknowns) by taking the maximum of corresponding entries
mergeLiftInfo = method();

mergeLiftInfo ( LiftInfo, LiftInfo ) := ( firstliftInfo , secondLiftInfo )->
(
    maxLiftDepth        := max ( firstliftInfo#"maxLiftDepth",       secondLiftInfo#"maxLiftDepth" );
    maxLatticeDimension := max ( firstliftInfo#"maxLatticeDimension",secondLiftInfo#"maxLatticeDimension");
    requiredLatticeDimension := null;  -- final requiredLatticeDimension != null only if it is present in both lift statistics!
    if ( firstliftInfo#"requiredLatticeDimension"=!=null and secondLiftInfo#"requiredLatticeDimension"=!=null ) then 
        requiredLatticeDimension =  max (firstliftInfo#"requiredLatticeDimension", secondLiftInfo#"requiredLatticeDimension");
    return createLiftInfo( maxLiftDepth, maxLatticeDimension, requiredLatticeDimension );
)


 
-- todo : introduce minimalLatticeDim and detect minimalLatticeDim by looking at LLLOutput in case a solution was found?

tryLatticeReduction = method (Options=>{"reductionOpts"=>new LiftOptions});



tryLatticeReduction (RingElement, Matrix, Matrix ) := LatticeReductionResult => opts-> (unknown, liftResult, nextLiftResult ) ->
(
    
    reductionOpts := opts#"reductionOpts";

    -- "TODO: check fails " assert(class reductionOpts === LiftOptions);
   

    foundMinPolynomialCandidate := false;
    latticeDimIncrementFkt := reductionOpts#"latticeDimIncrementFkt";
    verbose := reductionOpts#"verbose";
      --debug := reductionOpts#"debug";
    currentLatticeDim  := reductionOpts#"initialLatticeDim";
    minColumnNormDistanceFactor:= reductionOpts#"minColumnNormDistanceFactor";
    
    if verbose then print reductionOpts;
    latticeBasisVectorNormsList := {};
    lastColumnNormMin := -1;
    reducedLatticeBasis:= null;
    while (currentLatticeDim <= reductionOpts#"maxLatticeDim" ) do 
    (
        LLLInput := constructLLLInputFromLift(unknown, liftResult, currentLatticeDim );
        if (verbose) then 
        (print "time LLL:"; time reducedLatticeBasis = LLL( LLLInput, Strategy=>reductionOpts#"LLLStrategy" );) 
        else  reducedLatticeBasis = LLL( LLLInput , Strategy=>reductionOpts#"LLLStrategy");
        --
        -- test, if a solution have been found in this step (=foundMinPolynomialCandidate):
        bvec :=   sub(  matrix{apply(currentLatticeDim, i->unknown^i)}    ,    sub(nextLiftResult,ZZ) ) ;
        assert ( (ring bvec)===ZZ );

        columnNormData := computeNormalizedColumnNormsData(reducedLatticeBasis);

        if (verbose) then
        (
            print "columnNormData#max";
            print columnNormData.max;
            print "columnNormData#min";
            print columnNormData.min;
        );
        
        if ( sub(bvec*reducedLatticeBasis_{0}, ring nextLiftResult) == 0 and 
            (columnNormData.max//columnNormData.min)>minColumnNormDistanceFactor ) then 
            foundMinPolynomialCandidate = true;
        -- sometimes first condition (sub(bvec*reducedLatticeBasis_{0}, ring nextLiftResult) == 0) passes, but we do not have a solution. 
        -- Due to HC if we will use a higher lift, this could be detected at the end.
        ---
        -- stop condition for increasing lattice dimension: 
        --   for a generic stop condition it should be sufficient to get as input previous and current reducedLatticeBasis?
        
        latticeBasisVectorNormsList = latticeBasisVectorNormsList | { columnNormData } ;
        if (verbose) then  
            ( print "column norms: "; print ( toString columnNormData.normalizedNorms  | ", " | toString columnNormData.min); );
        
        if ( ( lastColumnNormMin==columnNormData.min ) or foundMinPolynomialCandidate ) then  
        (
            if (verbose) then  (print ("currentLatticeDimension: " |  currentLatticeDim););
            break;
        );
        lastColumnNormMin = columnNormData.min;
        assert( (latticeDimIncrementFkt ( currentLatticeDim )) > currentLatticeDim or currentLatticeDim == infinity );
        currentLatticeDim = latticeDimIncrementFkt(currentLatticeDim);
    );
    -- todo: Typ für Rückgabe einfuehren?
    result := new MutableHashTable;
    result.foundMinPolynomialCandidate= foundMinPolynomialCandidate;
    result.reducedLatticeBasis= reducedLatticeBasis;
    result.latticeBasisVectorNormsList= latticeBasisVectorNormsList;
    result.currentLatticeDim= currentLatticeDim;

    if (verbose) then 
    (
        print "tryLatticeReduction result:";
        print result;
    );
    return new LatticeReductionResult from result;
)


-- latticeBasisToPolynomial : use the first column 'latticeBasisVector' of the latticeBasis for a polymial construction:
--  return polynomial = latticeBasisVector_i*variable^i, i==0..#latticeBasisVector-1
latticeBasisVectorToPolynomial = method();
latticeBasisVectorToPolynomial(Matrix,RingElement) := RingElement => (latticeBasisVector,variable)->
(
    localVariable:=null;
    x:=null;
    if (variable===null) then 
    (
        localVariable = GlobalInternalPadicLiftResultVariable;
    )
    else
    (
        if (#select(gens ring variable, (curr)->(return curr===variable))==0) then
            error "variable is not a generator of a ring";
    
        localVariable =variable;    
    );
    nrows := numRows(latticeBasisVector);
    pol := matrix{apply(nrows , i->localVariable^i) }*(sub(latticeBasisVector,ring localVariable));
    return pol_0_0;
)


-- note : tryLatticeReduction wiped out from    computeMinPolys . Advances: Code readability. 
--        Potential drawbacks:  (implicit) copying   latticeBasis (return value) could issue memory  problems . 
-- Solution idea: return latticeBasis encapsulated in a HashMap  
-- computeMinPolys : 
-- tries to lift one of the unknown of a solution over a finite field for a system of equations( given as an ideal 'IZ' over ZZ) to characteristic zero 
--  using p-adic lifting and appying the Lenstra–Lenstra–Lovász(LLL)-algorithm to the lift.
--
-- parameter: equations ideal, solution, unknown (elements of gens ring coeffs or linear combinations) and self-explaining optional parameters
-- TODO: timeout as parameter...
-- TODO: Type safety?

computeSingleMinPolyEx = method (Options=>{"options"=>new LiftOptions});
computeSingleMinPolyEx (Ideal, Matrix, RingElement,RingElement) := 
ReducedPadicLiftResult => opts -> (IZ, solution, unknown, resultPolynomialVariable ) ->
(    
    --
    -- checking function parameters:
    checkLiftOptions (opts#"options");
    if not ( ring IZ === ring unknown) then error "an unknown is not element of the equations ideal ring";
    if not ( sub(IZ, solution) == 0)   then error "the given solution is not an element of the equations ideal vanishing set";

    liftOptions:= opts#"options";
    initialLiftDepth := liftOptions#"initialLiftDepth"; 
    maxLiftDepth    := liftOptions#"maxLiftDepth";
    verbose := liftOptions#"verbose";

    jacobianIZ := jacobian IZ;
    --
    currLiftDepth := 0;
    reducedLatticeBasis := null;
    liftResult := solution;
    latticeBasisVectorNormsList := null;
    nextLiftResult :=     nextLift( IZ, jacobianIZ, liftResult);
    foundMinPolynomialCandidate := false;
    currentLatticeDim := -1;
    --
    -- increase lift depth and perform LLL until a solution is found or maxLiftDepth is reached.
    while (currLiftDepth <= maxLiftDepth ) do (
        if (verbose) then print ( "--\n" | "currLiftDepth: " |  currLiftDepth);        
        --
        -- perform LLL onliy if (currLiftDepth >= initialLiftDepth ). 
                -- The condition is useful if minimalLiftDepth (=initialLiftDepth ) is known (e.g. from experimenting)
        if (currLiftDepth >= initialLiftDepth ) then    
        (
            lllReductionResult:= tryLatticeReduction(   unknown, liftResult,  nextLiftResult, "reductionOpts"=>liftOptions );
            
            foundMinPolynomialCandidate = lllReductionResult.foundMinPolynomialCandidate;
            currentLatticeDim =  lllReductionResult.currentLatticeDim;
            latticeBasisVectorNormsList = lllReductionResult.latticeBasisVectorNormsList;
            reducedLatticeBasis =  lllReductionResult.reducedLatticeBasis;
            if ( foundMinPolynomialCandidate  ) then  (
                --
                if (verbose) then print (" finalLiftDepth: " |  currLiftDepth);
                break;
            );
        );
        currLiftDepth = currLiftDepth+1;    

        liftResult = nextLiftResult;
        nextLiftResult =     nextLift(IZ, jacobianIZ, liftResult);
    );

    if (verbose) then print "reducedLatticeBasis";
    if (verbose) then print reducedLatticeBasis;

    result := new MutableHashTable;
    result.foundMinPolynomialCandidate = foundMinPolynomialCandidate;
    --result.reducedLatticeBasis = reducedLatticeBasis; --for debugging...
    --result.latticeBasisVectorNormsList= latticeBasisVectorNormsList; --for debugging...

    if ( not foundMinPolynomialCandidate  ) then
    (   
        result.polynomial = null;
        result.liftInfo = createLiftInfo( currLiftDepth, currentLatticeDim, null );
        result.unknown =unknown;
        return new ReducedPadicLiftResult from result;
    );   
    
    prepolynomial:=  latticeBasisVectorToPolynomial(reducedLatticeBasis_{0},resultPolynomialVariable);
    prepolynomialFactors := removeConstantFactors(factor prepolynomial);
    
  
    if ( #prepolynomialFactors>1) then
    {
       localLiftOptions := copy opts#"options";
       localLiftOptions#"initialLatticeDim" = localLiftOptions#"initialLatticeDim"-#prepolynomialFactors+1;
        return computeSingleMinPolyEx(IZ, solution, unknown, resultPolynomialVariable ,"options"=>localLiftOptions);
    };
    result.polynomial= value prepolynomialFactors;
    result.liftInfo = createLiftInfo( currLiftDepth, currentLatticeDim, (degree result.polynomial)#0 + 1 );
    result.unknown =unknown;
    return new ReducedPadicLiftResult from result;
)


-- ist für die Dokumentation doof. Führe optionen ein, mit dem Hinweis, dass (symbol option)=> benutzt werden soll um eine option zu setzten




testComputeSingleMinPolyEx=()->
(
    x:=null; x=symbol x;
    prime := 7;

    RQ := QQ[x];
    x = (gens(RQ))#0;
    FQ := (1/3*x+1/2) ; 
    IFQ := ideal FQ;

    -- reformulate with integer coeff ring to get correct results:
    IFZ := disposeRationalCoeffs(IFQ);  

    x = (gens ring IFZ )#0;
    solutionQQ := matrix{{-3/2}} ;
    solution := sub( solutionQQ ,ZZ/prime) ;
    result := computeSingleMinPolyEx( IFZ, solution, x,x);

    assert(result#polynomial==2*x+3);
)




-- todo: nur das Polynom oder auch die Gitterbasis-Matrizen (zu debugzwecken) zurueckgeben?

computeMinPolys = method (Options=>{"options"=>new LiftOptions});
computeMinPolys (Ideal, Matrix, List) := opts->(equationsIdeal, solutionPoint, unknownList)->
(
    checkLiftOptions (opts#"options");

    solutionChar := char ring solutionPoint;
    assert(solutionChar>0);

    SZ := ring equationsIdeal;
    --betti gens equationsIdeal

    varsnum := #unknownList;

    unknown :=null; 
    liftResult  :=null; 
    minimalPolynomialsTable := new MutableHashTable;


    liftInfo := createLiftInfo(0,0,0);


    apply(varsnum, currVarIdx->(
        unknown = (unknownList)#currVarIdx;

        -- 
        polResVar := GlobalInternalPadicLiftResultVariable;

        -- if the current unknown is a generator of the equations ideal ring, use the unknown as indeterminate for the minimal polynomial.
        -- otherwise the global GlobalInternalPadicLiftResultVariable is used. 
        if (#select(gens ring unknown, (curr)->(return (curr==unknown);))>0) then
            polResVar = unknown;

    
        -- adjust lift options.
        localOptions := opts#"options";
        
        if (localOptions#"initialLiftDepth"<liftInfo#"maxLiftDepth") then 
            localOptions#"initialLiftDepth"=liftInfo#"maxLiftDepth";

        if (localOptions#"initialLatticeDim"<liftInfo#"requiredLatticeDimension") then 
            localOptions#"initialLatticeDim"=liftInfo#"requiredLatticeDimension";
        
        -- todo: investigate: if a local variable 'localOptions' is defined, (symbol localOptions) has different behaviour
        liftResult = computeSingleMinPolyEx( equationsIdeal, solutionPoint, unknown, polResVar, "options"=>localOptions );
        liftInfo = mergeLiftInfo(liftInfo, liftResult.liftInfo );
        minimalPolynomialsTable#unknown = ( polResVar, liftResult.polynomial );
    ));
    
    return (minimalPolynomialsTable, liftInfo);
)


computeSingleMinPoly = method (Options=>{"options"=>new LiftOptions});
computeSingleMinPoly (Ideal, Matrix, RingElement) :=  ReducedPadicLiftResult =>  opts->(IZ, solution, unknown ) ->
(
    --return computeMinPolys(IZ, solution, unknown , "options"=>new LiftOptions); 

      polResVar := GlobalInternalPadicLiftResultVariable;

        -- if the current unknown is a generator of the equations ideal ring, use the unknown as indeterminate for the minimal polynomial.
        -- otherwise the global GlobalInternalPadicLiftResultVariable is used. 
        if (#select(gens ring unknown, (curr)->(return (curr==unknown);))>0) then
            polResVar = unknown;

    return computeSingleMinPolyEx(IZ, solution, unknown , polResVar, "options"=>opts#"options");
)


doc ///
    Key
        computeSingleMinPoly
        (computeSingleMinPoly ,Ideal,Matrix,RingElement)
    Headline
        compute a lift to a extension  of  QQ for a single unknown 
    Inputs
        equationsIdeal: Ideal 
            the equations ideal ( only integer coefficient ring is supported )
        solution: Matrix
            an element of the ideal vanishing set over a prime field
        unknown: RingElement
            the unknown to lift (element of 'gens ring {\tt equationsIdeal}')
--        resultVariable: RingElement
--            the result is a polynomial in {\tt resultVariable}; the roots of the polynomial are solutions for {\tt unknown}.
    Outputs
        : ReducedPadicLiftResult
            the minimal polynomial data with some computation statistics for the unknown.
    Usage
         liftResult = computeSingleMinPoly(equationsIdeal,solution,unknown)
         liftResult = computeSingleMinPoly(equationsIdeal,solution,unknown,"options"=>liftOptions)
    Description
        Text
            Given an equations ideal with integer coefficients and an element {\tt solutionPoint } of the vanishing set over a finite field, {\tt computeMinPolys} tries to compute the minimal polynomial for a given unknown. \break  \break
            The function contains a loop of two steps: p-adic Lift and LLL reduction. \break
            The process continues until the computation succeeded or stopping conditions are reached,
            see @TO LiftOptions@.       
        Example          
        Text
           \break  Example: find (subset of the) vanishing set for an ideal (here IFZ)
        Example      
            RQ = QQ[x,y];
            FQ = 33/4*x^3+19/4*x^2-81/4*x-1;          
            IFQ = ideal FQ    
            IFZ = disposeRationalCoeffs(IFQ)
            RZ  = ring IFZ;
        Text
           \break the solutions over a finite field can be found via "brute force" -  omitted here
        Example          
            K = ZZ/11;
            solution1 = matrix{ {1_K , 1_K} }; solution2 = matrix{ {5_K, 1_K} };
        Text
           \break check solutions
        Example
            assert( sub(IFZ, solution1) ==0 );   assert( sub(IFZ, solution2 ) ==0);
        Text        
           \break check if padic lifting is possible at all : jacobian has to be invertible.
        Example          
            JFZ = jacobian IFZ;
            JFp = sub(JFZ,solution2);
            assert(rank JFp==numColumns JFp);    
        Text    
            \break do it
        Example          
            liftResult2 = computeSingleMinPoly (IFZ, solution2, (gens ring IFZ)#0 );
        Text
            \break check lifted result correctness
        Example          
            solutionIdeal2 = ideal substitute(liftResult2.polynomial,(gens ring liftResult2.polynomial)#0=>(gens RZ)#0);
            assert ( isSubset(IFZ,solutionIdeal2) );                 
        Text
           \break the roots of the liftResult2.polynomial are solutions for FZ=0 , \break
            in this example this can be seen by looking at factored FZ
        Example          
            liftResult2.polynomial
            factor FQ   
    Caveat
        as coefficientring currently only ZZ is supported and the solution set has to be 0-dimensional. 
///

testComputeMinPolys = ()->
(
    x := null; x = symbol x;
    y := null; y = symbol y;
    prime := 7;

    RQ := QQ[ x,y ];
    x = (gens(RQ))#0;
    FQ := { (1/3*x+1/2), ((y-1/2) ) } ; 
    IFQ := ideal FQ;
    -- reformulate with integer coeff ring to get correct results:
    IFZ := disposeRationalCoeffs(IFQ);  

    x = (gens ring IFZ )#0;
    y = (gens ring IFZ )#1;
    solutionQQ := matrix{{-3/2,1/2}} ;
    solution := sub( solutionQQ ,ZZ/prime) ;
    result := computeMinPolys( IFZ, solution, { x, y } );

    assert(result#0#x#1==2*x+3);
    assert(result#0#y#1==2*y-1);

    result = computeMinPolys( IFZ, solution, { y, x  } );

    assert(result#0#x#1==2*x+3);
    assert(result#0#y#1==2*y-1);

    result = computeMinPolys( IFZ, solution, {  x  } );
    assert(result#0#x#1==2*x+3);

    -- a test, where startingLatticeDim is to big
    liftOptions := new LiftOptions;
    liftOptions#"initialLatticeDim"=6;
    result = computeMinPolys( IFZ, solution, {  x  },"options"=>liftOptions );
    assert(result#0#x#1==2*x+3);

    liftOptions  = new LiftOptions;
    liftOptions#"initialLatticeDim"=3;
    result = computeMinPolys( IFZ, solution, {  x  },"options"=>liftOptions );
    assert(result#0#x#1==2*x+3);


    result = computeMinPolys( IFZ, solution, {    } );

    -- root gluing test:
    result = approxComplexSolutions( IFZ, solution  );
)




--computeMinPolys (Ideal, Matrix, List) := (equationsIdeal, solutionPoint, unknownList)->
--(
--    return computeMinPolys(equationsIdeal, solutionPoint, unknownList, new LiftOptions)
--)

doc ///
    Key
        computeMinPolys
        (computeMinPolys ,Ideal,Matrix,List)
    Headline
        compute a lift to (an extension field of) rationals
    Inputs
        equationsIdeal: Ideal 
            the equations ideal ( only integer coefficient rings supported )
        solutionPoint: Matrix
            an element of the {\tt equationsIdeal} vanishing set over a prime field
        unknownList: List
            the unknowns  of interest (element of {\tt gens ring equationsIdeal} )
        "options"=> LiftOptions
            . If computation takes too long, consider to customize {\tt options }.
    Outputs
        : Sequence
            a pair   { \tt  (polynomialHashTable, }   @TO LiftInfo@ { \tt ) } . 
            Assume 'unknown' is an element of the input parameter {\tt unknownList}. \break
            Then roots of polynomial {\tt polynomialHashTable#unknown#1} ( with polynomial variable {\tt polynomialHashTable#unknown#0} )
            are solutions for the 'unknown'  .
    Usage
         ( polynomialHashTable, liftInfo ) = computeMinPolys( equationsIdeal, solutionPoint, unknownList);
         ( polynomialHashTable, liftInfo ) = computeMinPolys( equationsIdeal, solutionPoint, unknownList, "options"=>options);
    Description
        Text
            Given an equations ideal with integer coefficients and an element {\tt solutionPoint } of the vanishing set over a finite field, 
            {\tt computeMinPolys} tries to lift the {\tt solutionPoint } to (an extension field of) rationals. \break  \break
            The function contains a loop of two steps: p-adic lift and lattice reduction. \break
            The process continues until the computation succeeds or  stopping conditions (default or custom) are reached,
            see @TO LiftOptions@. \break
            The  roots of returned polynomials are solutions for the corresponding unknows.\break\break
            See also @TO computeSingleMinPoly@.
        Example          
        Text
           \break  Example: find (subset of the) vanishing set for an ideal {\tt IFZ }
        Example          
            RQ = QQ[symbol x,symbol y];
            FQ = (33*x^3+19*x^2-81*x-4)*y;          
            IFQ = ideal FQ;           
            IFZ = disposeRationalCoeffs(IFQ)
            RZ  = ring IFZ;
        Text
           \break Look at the problem over a finite field
        Example          
            K = ZZ/3; RK =K[symbol x,symbol y]; 
            IFK =sub(IFZ,RK)
         
        Text
           \break The solutions over a finite field can be found via "brute force": \break
        Example
            allPoints = flatten apply((char RK),x->apply(char RK, y->matrix {{ x_K, y_K }})); 
            --select the vanishing set of ideal IFK from all points in RK -- stand vorher als Text; finde ich aber unnoetig.
            solutions = select (allPoints,(point)->( sub( IFK,point) == 0) )  
        Text        
             \break Select points capable for lifting : jacobian has to be at least pseudoinvertible.
        Example 
              liftCapablePointList = select (solutions,(point)->( JFp = sub( jacobian IFZ,point); rank JFp==numColumns JFp )) 
        Text        
             \break Select a point for lifting - all in {\tt liftCapablePointList} are suitable!
        Example 
              solutionPoint = liftCapablePointList#1
--        Text    
--            \break Select liftable unknown candidates : jacobian has to be at least pseudoinvertible.
--        Example          
--            capableUnknownList = select(gens ring IFZ, (unknown)->( sub (diff(unknown, FZ), solutionPoint)!=0 ) )
--            unknown = capableUnknownList#0;
        Text    
            \break Compute minimal polynomial for x:
        Example        
            unknown = (gens RZ)#0
            ( minimalPolynomialsTable , liftInfo )  = computeMinPolys  ( IFZ, solutionPoint, {unknown} );
        Text
            \break The roots of the polynomial {\tt minimalPolynomialsTable#unknown }  are solutions for variable {\tt unknown}. \break
            Compare with factored initial polynomial FZ!
        Example          
            minimalPolynomialsTable#unknown
            factor FQ  
        Text
            \break Optiional: check lifted result correctness without looking at the solution \break
            (whether {\tt IFZ} is contained in { \tt equationsIdeal } )
        Example          
            equationsIdeal = ideal  minimalPolynomialsTable#unknown#1;
            --equationsIdeal2 = ideal substitute(minimalPolynomialsTable#unknownId,(gens ring polynomials#unknownId)#unknownId=>unknownList#unknownId);
            isSubset(IFZ,equationsIdeal)
            assert ( isSubset(IFZ,equationsIdeal) );           
    Caveat
        As coefficient ring currently only integers (@TO ZZ@) are supported and the solution of interest has to be isolated and smooth. \break
        It is possible that not the smallest appropriate latticeDimension will be used for reduction.
///


approxComplexSolutions = method (Options=>{"options"=>new LiftOptions,"decimalPrecision"=>10});
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
idealPointsApproxData=( inputIdeal, solutionPoint, minPolyTable, approxSolutions, mergedLiftInfo, unknownList )->
(
    approxSolutionData := new MutableHashTable;
    approxSolutionData#"inputIdeal"  =  inputIdeal ;
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



-- may have problems in case unknownList is different from (gens ring inputIdeal order). Therefore removed the parameter for a moment.
approxComplexSolutions (Ideal, Matrix) := opts->(inputIdeal, solutionPoint)->
(   

    unknownList := gens ring inputIdeal;

    if (not pariGpIsPresent() ) then 
        error "please install Pari/GP or build Macaulay2 with pari (add 'pari' to the --enable-build-libraries configure parameter. ) to use 'approxComplexSolutions'";

    --localDecimalPrecision := opts#"options"#"decimalPrecision";
    localDecimalPrecision := opts#"decimalPrecision";

    rootCalculator := opts#"options"#"rootCalculator";

    characteristic:=char ring solutionPoint;


    (minimalPolynomialsTable, mergedLiftInfo ) := computeMinPolys(inputIdeal, solutionPoint, unknownList, 
                                    "options" => opts#"options" );

    pairedRootRootList := new MutableList;

    opts#"options"#"initialLiftDepth"   = mergedLiftInfo#"maxLiftDepth"+1 ; --# is a heuristic. could be suboptimal for generic problems.
    opts#"options"#"initialLatticeDim"  = mergedLiftInfo#"requiredLatticeDimension" ;

    -- well, what should we do, return some info or return null in case something failed?
    if (mergedLiftInfo#"requiredLatticeDimension"===null) then  
    (
        return null;
        --return idealPointsApproxData( inputIdeal, solutionPoint, minimalPolynomialsTable, {}  , mergedLiftInfo, unknownList );
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
           (compatibilityMinimalPolynomials, compatLiftInfo ) := computeMinPolys(inputIdeal, solutionPoint, {newUnknown}, 
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
    
    return idealPointsApproxData(  inputIdeal, solutionPoint, minimalPolynomialsTable, preApproxSolutions  , mergedLiftInfo, unknownList );
    
    --return approxComplexSolutionData;
)


doc ///
    Key
        approxComplexSolutions
        (approxComplexSolutions ,Ideal, Matrix)
    Headline
        given a solution over a prime field for a system of equations, compute the corresponding complex solutions
    Inputs
        equationsIdeal: Ideal 
            the equations ideal ( only integer coefficient ring is supported )
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
         complexSolutionsData = approxComplexSolutions( equationsIdeal, solutionModPrime, unknowns )
         complexSolutionsData = approxComplexSolutions( equationsIdeal, solutionModPrime, unknowns,"options"=>options, "decimalPrecision"=> decimalPrecision )

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
            equationsIdeal =  disposeRationalCoeffs IFQ
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
    IFZ := disposeRationalCoeffs(IFQ);  

    x = (gens ring IFZ )#0;
    y = (gens ring IFZ )#1;
    solutionQQ := matrix{{-3/2,1/2}} ;
    solution := sub( solutionQQ ,ZZ/prime) ;


    -- root gluing test:
    result := approxComplexSolutions( IFZ, solution  );
)
-------------


------------------------------------------------------------------------------------------------------------
--Example from finite field experiments (Bothmer, et al)
finiteFieldExperimentsExampleTest =()->
(
    x := null; x = symbol x;
    y := null; y = symbol y;
    R := ZZ[x,y];
    x = (gens R)#0;
    y = (gens R)#1;
    -- the equations
    I := ideal (-8*x^2-x*y-7*y^2+5238*x-11582*y-7696,
    4*x*y-10*y^2-2313*x-16372*y-6462);

    prime := 7;
    Fp := ZZ/prime;

    P1 := matrix { { 2_Fp, 3_Fp } };
    P2 := matrix { { 5_Fp, 5_Fp } };

    assert( sub(I,P1)==0 );
    assert( sub(I,P2)==0 );    

    liftedData := computeMinPolys(I,P1,{x,y});
    assert( liftedData#0#x#1== x-1234 );
    assert( liftedData#0#y#1== y+774 );
    computeMinPolys(I,P2,{x,y});
)


TEST ///
debug padicLift
padicLiftProtect()
testLiftStep()
///


TEST ///
debug padicLift
padicLiftProtect()
testLiftPoint()
///

TEST /// 
debug padicLift
padicLiftProtect()
testRemoveConstantFactors()
///

TEST ///
debug padicLift
padicLiftProtect()
testDisposeRationalCoeffs()
///


TEST ///
debug padicLift
padicLiftProtect()
testNestedRingCoeffsLCMDenominator()
///
         

TEST ///
debug padicLift
padicLiftProtect()
testTensoredDisposeRationalCoeffs()
///

TEST ///
debug padicLift
padicLiftProtect()
testComputeSingleMinPolyEx()
///

TEST ///
debug padicLift
padicLiftProtect()
testComputeMinPolys()
///

TEST ///
debug padicLift
padicLiftProtect()
finiteFieldExperimentsExampleTest()
///


TEST ///
debug padicLift
padicLiftProtect()

if pariGpIsPresent() then 
  testComputeRootsWithGP();

///

TEST ///
debug padicLift
padicLiftProtect()
  testExtractArrayString();

///



-- needs a test where gluing is tested! - done
-- needs a test, where startingLatticeDim is to big - done
-- needs a test, where the 
-- needs tests, where nonstandard gluing cases, as described in our paper are tested. 


----------------------------------------------------- sandbox ---------------------------------------------------

--padicLiftInfo =  ( infoLevel, message)=>
--(
--    if ( padicLiftVerboseLevel >= infoLevel ) then
--    print("--padicLiftInfo: "| toString(message)|"\n");
--)


testOptions = method (Options=>{ protect optionName=>false})
testOptions (ZZ)  := ZZ => opts1 ->  (param) ->
(
    print opts1;
    print opts1#optionName;
    return 1_ZZ;
)


-- playing with testing the (ugly) gp root computation routines of this package
prepareCompareCommand=(initialNumString,convertedNumString,decimals)->
(
    command:= "\\p " | decimals |";\n";
    command=command |"print (\"beginOutput\");\n";
    command =command |"\n if (" |initialNumString |"==" | convertedNumString |", print (\"true\"))\n";
    command =command |"\n if (" |initialNumString |"!=" | convertedNumString |", print (\"false\"))\n";
    return command;
)

numbersAreEqual = (initialNumString,convertedNumString,decimals)->
(
    command := prepareCompareCommand(initialNumString,convertedNumString,decimals);
    resultStr := callGP(command);
    pos1 :=regex("beginOutput",resultStr);
    assert(pos1=!=null);
    subs := substring(pos1#0#0,resultStr);
   -- subs3 = substring(0,3,subs);
    pos1=regex("false",subs);
    if (pos1=!=null) then 

        return false;
    pos1=regex("true",subs);
    if (pos1=!=null) then 
        return true;
    assert(false);
)

------------------------------------------------------------
computeMinPolys
