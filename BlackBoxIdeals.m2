newPackage(
     "BlackBoxIdeals",
     Version => "0.2", 
     Date => "23.02.2017",
     Authors => {{
           Name => "Jakob Kroeker", 
           Email => "jakobkroeker.academic@spaceship-earth.net", 
           HomePage => "http://www.crcg.de/wiki/User:Kroeker"},{
           Name => "Hans-Christian Graf v. Bothmer", 
           Email => "bothmer@math.uni-hannover.de", 
           HomePage => "http://www.crcg.de/wiki/Bothmer"}    
      },
     Configuration => {},
     PackageExports => {"M2Logging"},
     Headline => "black boxes for explicit and implicitly given ideals",
     DebuggingMode => true,
     CacheExampleOutput => false,
     AuxiliaryFiles=>true
)


needsPackage "M2Logging";



export { 
    "OnComponentAnswerStrategy",
    "PlainTextSmoothnessInfoWithAnswerPair",
    "NullIfNotSmoothStrategy",
    "ExceptionIfNotSmooth",
    "SmoothnessInfoWithAnswerPair",
    "plainTextSmoothnessInfoWithAnswerPair",
    "nullIfNotSmoothStrategy",
    "exceptionIfNotSmooth",
    "smoothnessInfoWithAnswerPair",
    "refineInterpolation",
    "maximalConditions",
    "monomialBasisSize",
    "basicJetLengthHeuristic",
    "constantJetLengthHeuristic",
    "setBlackBoxLogLevel",
    --"compatible",
    "continueJet",
    "continueJetOrException",
    "JetLengthHeuristic",
    "ConstantJetLengthHeuristic",
    "BasicJetLengthHeuristic",
    "InterpolationMonomialDegreeHeuristic",
    "BasicInterpolationMonomialDegreeHeuristic",
    "ConstantInterpolationMonomialDegreeHeuristic",
    "sameComponentAt",     
    "componentNamesAt",
    "bare",
    "jetStatsAt",
    --"locus",
    --"joinJetSets",
    "firstElem",
    "JetSet",
    "addElement",
    "isDerivedFrom",
    "Jet",
    "pointProperties",
    "memberMethods",
    "attributes",
    "clearCoeffDenominators",
    "BlackBoxIdeal",
    "BlackBoxParameterSpace",
    "blackBoxParameterSpace",
    "blackBoxIdeal",
    "blackBoxIdealFromEvaluation",
    "BlackBoxLogger",
    "getEpsRing",
    --"bestJetAt",
    "jetAt",
    "jetAtOrException",
    "jetAtWithInfo",
    "jetAtOrNull",
    "isCertainlySingularAt",
    "isProbablySmoothAt",
    "keysWithoutSymbols",    
    "guessAcceptedParameterNumber",
    "dropDegreeInfo",
    "createInterpolatedIdeal",
    "interpolateAt",
    "interpolate",
    "MapHelper",
    "InterpolatedIdeal",
    "SingularPointException",
    "PointNotOnBlackBox",
    "createMapHelper"
}


userDictHasKey = (key)->
(
    keyStr := toString key;
    userKeyList := apply(keys User#"private dictionary", b->toString b);
    pos := position(userKeyList, (v)->v==keyStr);
    return (pos =!= null);
);
 
userDictKey = (key)->
(
    keyStr := toString key;
    userKeyList := apply(keys User#"private dictionary", b->toString b);
    pos := position(userKeyList, (v)->v==keyStr);
    return (keys User#"private dictionary")#pos;
);


idealBlackBoxesProtect = ()->
(
    --protect eps;
    protect withChecks;
    protect enableChecks;
    protect disableChecks;
    protect interpolateComponents;
    protect jacobianAt;
    protect rankJacobianAt;
    protect valuesAt;
    protect unknownIsValid;
    protect numVariables;
    protect numGenerators;
    protect isZeroAt;
    protect registerPointProperty;
    protect rpp;
    protect numTrials;
    protect setSingularityTestOptions;
    protect singularityTestOptions;
    protect updateSingularityTest;
    protect setPointProperty;
    protect setValuesAt;
    protect checkInputPoint;
    protect deduceNumGenerators;
    protect setIsZeroAt;
  
     
    protect unknowns;
    protect pointProperty;
    protect updatePointProperty;
    protect setJacobianAt;
    protect equations;

    protect hasPointProperty;
    protect pointPropertiesAsSymbols;
    --protect memberMethods;
    --protect attributes;
    --protect knownProperties;
    --protect createMapHelper;
    protect updateBlackBox;
    protect setInterpolator;
);

--todo: fix dublicate code,  -  padicLiftProtect and padicLiftExport

idealBlackBoxesExport = ()->
(
    exportMutable("transformedAnswer");
    exportMutable("onComponentAnswerStrategies");
    exportMutable("setOnComponentAnswerStrategy");
    exportMutable("onComponentAnswerStrategy");
    exportMutable("onComponentPrecision");
    exportMutable("interpolateComponentAt");
    exportMutable("interpolateComponentsAt");
    exportMutable("interpolatedComponentsAt");
    exportMutable("interpolatedComponents");
    exportMutable("interpolatedComponentNames");
    exportMutable("interpolatedComponentByName");        
    --exportMutable("reset");
    exportMutable("setOnComponentPrecision");
    exportMutable("resetInterpolation");
    exportMutable("setInterpolator");
    exportMutable("setMonomialDegreeHeristic");
    exportMutable("componentsAt");
    exportMutable("interpolationTargetJetLength");
    exportMutable("setJetLength");
    exportMutable("setAdditionalJetLength");
    exportMutable("setJetLengthHeuristic");
    exportMutable("setMonomialDegreeHeuristic");
    exportMutable("targetMonomialDegree");
    exportMutable("componentNames");
    exportMutable("clearCache");
    exportMutable("minComponentDegree");
    exportMutable("maxInterpolationDegree");
    exportMutable("maxComponentDegree");
    exportMutable("blackBox");
    exportMutable("dropComponent");
    exportMutable("increaseInterpolationJetLength");
    exportMutable("decreaseInterpolationJetLength");
    exportMutable("setSameComponentPrecision");
    exportMutable("setComponentNamePrefix");
    exportMutable("sameComponentTargetJetLength");
    exportMutable("additionalJetLength");
    exportMutable("interpolator");
    exportMutable("componentByName"),
    exportMutable("renameComponent"),
    exportMutable("componentNameInUse"),
    exportMutable("componentNamesInUse"),
    exportMutable("jetSet");
    exportMutable("setName");
    --exportMutable("name");
    --exportMutable("name");
    exportMutable("isOnComponent");
    exportMutable("isOnInterpolatedComponent");
    exportMutable("enableChecks");
    exportMutable("disableChecks");
    exportMutable("withChecks");  
   
    exportMutable("interpolateComponents");
 
    exportMutable("eps");
    exportMutable("jacobianAt");
    exportMutable("rankJacobianAt");
        
    exportMutable("valuesAt");
        
    exportMutable("unknownIsValid");                 
    exportMutable("numVariables");  
    exportMutable("numGenerators");
    exportMutable("isZeroAt");      
    exportMutable("registerPointProperty"); 
    exportMutable("rpp");
    exportMutable("upp");
    exportMutable("setPointProperty");
    exportMutable("setValuesAt");    
    exportMutable("checkInputPoint");
    exportMutable("deduceNumGenerators");
    exportMutable("setIsZeroAt");
    
    exportMutable("unknowns");
    exportMutable("pointProperty");
    exportMutable("updatePointProperty");
    exportMutable("setJacobianAt");
    exportMutable("equations");

    exportMutable("hasPointProperty");
    exportMutable("pointPropertiesAsSymbols");
    --exportMutable("memberMethods");
    --exportMutable("attributes");
    exportMutable("numTrials");
    exportMutable("setSingularityTestOptions");
    exportMutable("updateSingularityTest");
    exportMutable("singularityTestOptions");
    exportMutable("updateBlackBox");
    
)


--load "./BlackBoxIdeals/Exceptions.m2";
------------------------------------------------
-- EXCEPTIONS
------------------------------------------------
SingularPointException = new Type of HashTable;


singularPointException = ()->
(
    return new SingularPointException from {}
)

PointNotOnBlackBox = new Type of HashTable;

------------------------------------------------
-- END EXCEPTIONS
------------------------------------------------
--load "./BlackBoxIdeals/Utils.m2";

------------------------------------------------
-- UTILS
------------------------------------------------
--
-- guessAcceptedParameterNumber(): 
--
-- find out for a function, how many parameters it does accept.
-- if a function accepts variable number of parameters, returns null
-- if it did not find 'numparms:' in the disasseble string, returns null.
--

guessAcceptedParameterNumber = method();

guessAcceptedParameterNumber( Function ) := ZZ => (foo)->
(
    lst := disassemble foo;

    --bblog.debug  ("disassemble result: " | lst );
    lst = separate( " ", lst );

    restargsPos := position( lst, (str)-> str=="restargs:" );
    numparmsPos := position( lst, (str)-> str=="numparms:" );

    if restargsPos=!=null then 
    (
        newlst := drop (lst, restargsPos);
        restargsPos = position( newlst, (str)-> str=="restargs:" ); 
        
        if restargsPos=!=null then 
        (
            --bblog.info ("do not know how to handle methods with a chain of several '->' ");
            --return null; 
        );
    );

    if numparmsPos===null then 
    (
        --bblog.warning (" warning: did not find the position of 'numparms:' in dissasseble string ");
        return null; 
    )
    else
        return  value lst#(numparmsPos+1);
);


-- guessAcceptedParameterNumber():
--
--   find out for a method, how many parameters it does accept.
--   Works only if a single function is installed for that method;
--   if multiple functions are installed for the same method, returns null.
--   status: beta.
--
guessAcceptedParameterNumber( MethodFunction ) := ZZ=>(foo)->
(
    func := apply( methods foo , m-> lookup m);
    if #func==1 then 
    (
        return  (guessAcceptedParameterNumber func#0);
    );
    if #func>1 then 
    (
        --bblog.info ("did not expect a method with multiple installed functions. ");
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
    -- (foo(isPrime))(4,3)

    foo = a->(a,b);

    assert( 1==guessAcceptedParameterNumber foo );

)


TEST ///
    debug BlackBoxIdeals
    idealBlackBoxesProtect()
    testGuessAcceptedParameterNumber()
///


-- polynomialLCMDenominator()
--
-- computes the least common multiple of the denominators of the polynomial (rational) cofficients;
-- that means if we have a polynomial = sum { (a_i/b_i)*monomial_i }, a_i and b_i integers,
-- then the function returns LCM( (b_i) )
--
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
        ( 
            error("expected rationals as coefficient ring!"); 
        );
    );
    if (coeffRng===QQ) then   LCMDenominator =  lcm apply(summands ,j-> denominator j ) ;
    return LCMDenominator;
)


-- clearCoeffDenominators() 
--
-- converts an ideal with rational coefficients 
-- to an ideal with integer coefficients while preserving the vanishing set.
-- e.g. if sub(IdealWithRationalCoeffs,point)==0, then  
-- sub( clearCoeffDenominators(IdealWithRationalCoeffs),point)==0 and vice versa
--
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



-- testClearCoeffDenominators()        
--
-- test conversion of an ideal with rational coefficients to an ideal 
-- with integer coefficients while preserving the vanishing set. 
--
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

-- testNestedRingCoeffsLCMDenominator()
-- 
-- test polynomialLCMDenominator ( computing the least common multiple 
-- of the denominators of polynomials with rational coefficients)
-- in case that the rationals(QQ) are not the coefficient ring
--
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

-- testTensoredClearCoeffDenominators()
-- 
-- test conversion of an ideal with rational coefficients to an ideal 
-- with integer coefficients while preserving the vanishing set. 
-- (clearCoeffDenominators)
-- the special case, that the ring the equations belong to is a tensor product of other rings. 
--
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
------------------------------------------------
-- END UTILS
------------------------------------------------


Jet = new Type of HashTable;

jetObject := method();

-- jetObject.parent is currently black box

jetObject( Thing, Thing, Thing, ZZ) := Jet=>(parent, point, value, jetLength)->
(
      jetObj := new HashTable from { ("value",value),
                                     ("jetLength",jetLength),                                   
                                     ("parent",parent),
                                     ("point",point),
                                   };
    jetObj = newClass (Jet, jetObj);
    return jetObj;
);

 

net (Jet) := Net =>(jet)->
(
    sss :=  { --( net( "parent")  | " <= " | net ( class jet#"parent") ),
              --( net( "point")  | " <= " | net (jet#"point") ), 
              ( "(" | net (jet#"jetLength") |") " | net (jet#"value") )
              --( net( "length")  | " = " | net (jet#"jetLength") )              
            } ;    
    return net stack sss;
)

-- TODO naming: length, or precision?
length (Jet) := ZZ =>(jet)->
(
    return jet#"jetLength";
);

TEST /// 
    -- test 'length(Jet)'
    Fp = ZZ/101
    R = Fp[x,y]
    I = ideal(x^2-y^2+x^3)
    bbI = blackBoxIdeal I;
    point = matrix{{3,6_Fp}}
    bbI.isZeroAt(point)
    j = jetAt(bbI,point,3)
    assert(length j == 3)
///


sub (Ideal, Jet ) := Thing =>(I,jet)->
(
    return (sub(I,jet#"value"));
)



-- swith between protect and export - both are not possible!

--idealBlackBoxesProtect() -- protect the symbols for checking package correctness: no symbol variables should be overwritten by accident!
idealBlackBoxesExport(); -- export the symbols to make the package work 


needsPackage "SimpleDoc";
needsPackage "Text";


BlackBoxLogger = Logger("BlackBoxIdeals");

-- todo: how to switch this on and off by the user? using closures again? 
--if BlackBoxIdeals#Options#DebuggingMode then 
--    BlackBoxLogger.setLogLevel(LogLevel.DEBUG);

if BlackBoxIdeals#Options#DebuggingMode then
    errorDepth=0;



bblog := BlackBoxLogger; --shortcut

--if BlackBoxIdeals#Options#DebuggingMode then bblog.setLogLevel(LogLevel.DEBUG);

setBlackBoxLogLevel = (level)->
(
    BlackBoxLogger.setLogLevel(level);
);




savedEpsRings := new MutableHashTable;


-- getPropertySymbols() 
-- 
-- Returns all relevant symbols (propertySymbols) 
-- from different M2 scopes corresponding to a given method name (propertyName).
-- Background: suppose the user registered a point property with name 'PP' 
-- in a black box 'bb' and want to access it using the  '.' operator.
-- To achieve that, internally we have to add in the HashTable 'bb' an entry with 'symbol(PP)' as a key 
-- and the point property function as a value.
-- Now in Macaulay2  symbols are bound to scopes and thus in different scopes symbols with same name may coexist.
-- So it may happen that the registered key symbol 'PP' differes from the to the user visible symbol 'PP' at the top level
-- and the access 'bb.PP' leads to 'error: key not found in hash table'
-- The workaround is to register the same property at least twice using the symbol  getGlobalSymbol( 'PP' ) as a key.
-- I'm not sure, but for some situations it was also necessary to use the symbol from 'User#"private dictionary"'
-- via getGlobalSymbol( User#"private dictionary",'PP')
-- and the symbol from the package dictionary "BlackBoxIdeals.Dictionary"
-- see also https://github.com/jakobkroeker/FiniteFieldExperiments.M2/issues/116
-- 
-- TODO question: what happens in case the user loads a package which defines the same symbol?

-- artificial example:  
-- restart
-- loadPackage "BlackBoxIdeals"
-- fragileHT = new HashTable from {
--   getGlobalSymbol(User#"private dictionary", "valuesAt") => 5
-- }
-- fragileHT.valuesAt -- error: key not found in hash table
-- 

getPropertySymbols := method ();

getPropertySymbols(String) := List => (propertyName)->
(
    propertySymbols := {} ;
    try (  propertySymbols = propertySymbols | { getGlobalSymbol(BlackBoxIdeals.Dictionary, propertyName)} );

    -- todo question: should the symbol in the users private dictionary always be created?
    try  (  propertySymbols = propertySymbols | { getGlobalSymbol propertyName} ) else 
    ( 
        propertySymbols =  propertySymbols | { getGlobalSymbol( User#"private dictionary", propertyName); }
    );
    return propertySymbols;
);


--
-- package-global symbol for "eps"
--
--geps := getSymbol "eps";

geps := symbol eps;

isDerivedFrom = method();

isDerivedFrom (Thing, Type) := Boolean => (testee, expectedType)->
(
    if (class testee===expectedType) then return true;
    
    parentClass := parent class testee;
    
    while (parentClass =!= Thing) do
    (
        if (parentClass===expectedType) then return true;
        parentClass = parent parentClass;
    );
    if (parentClass===expectedType) then return true;    
    return false;
)

getEpsRingFunction := (coeffring, epsDim)->
(
    assert (isDerivedFrom(coeffring,Ring));
    assert (isDerivedFrom(epsDim,ZZ) or isDerivedFrom(epsDim,InfiniteNumber));
    
    if isDerivedFrom(epsDim,InfiniteNumber) then
    (
       if epsDim=!=infinity then error ("epsDim " | toString epsDim |"not supported");
    );
    if isDerivedFrom(epsDim,ZZ) then
    (
        if (epsDim<0) then error("expected epsDim>0 ");
    );
    
    leps := geps;   
    epsRng := null;
    eps := null;
    
    if not (savedEpsRings#?(coeffring,epsDim) ) then 
    (
        polRing := coeffring[leps];
        leps = first gens(polRing);
        if epsDim===infinity then 
        (
            savedEpsRings#(coeffring,epsDim) = polRing            
        )
        else
        (
            savedEpsRings#(coeffring,epsDim) = polRing/(leps^(epsDim+1));            
        );
        epsRng = savedEpsRings#(coeffring, epsDim);
        epsRng#"epsPrecision" = epsDim;
        eps = first gens epsRng;
        (savedEpsRings#(coeffring,epsDim))#"eps" = eps;
        (savedEpsRings#(coeffring,epsDim)).eps = eps;        
        --for symb in getPropertySymbols("eps") do 
        -- (
        --   assert(symb=!=null);
        --   (savedEpsRings#(coeffring,epsDim))#symb  = eps;
        --)
    ); 
    epsRng = savedEpsRings#(coeffring, epsDim);
    eps = (gens epsRng)#0;
    return epsRng
)



getEpsRing = method();

-- todo: get rid of duplicate code in getEpsRing()...

getEpsRing(Ring, InfiniteNumber) := Ring => (coeffring, epsDim)->
(
    return getEpsRingFunction(coeffring, epsDim);   
)


getEpsRing(Ring, ZZ) := Ring => (coeffring, epsDim)->
(
    return getEpsRingFunction(coeffring, epsDim);
)


testEpsRing = ()->
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



-- deduceNumGenerators():
--
--   for an (ideal) blackbox, determine the number of   generators (or equations). This is possible in case 
--   the blackbox provides the 'valuesAt'-property.
--
--   todo: this is a generic version. if the blackbox is given by an ideal I,
--        the generators can be determined easily by #(gens ideal I)
--   todo: what should the user do, if deduceNumGenerators() fails? 
--
deduceNumGenerators := (blackBox)->
(
    bblog.debug(" enter deduceNumGenerators ") ;

    if not blackBox.hasPointProperty("valuesAt") then 
    ( 
        error ("cannot determine the number of generators/equations, since valuesAt is missing");
    );
    try ( numVar:=blackBox.numVariables; ) else (  error (" blackBox.numVariables is missing"); );
    
    computed  := false; 
    maxTrials := 100;
    currTrial := 0;
    rng := blackBox.coefficientRing;
    numGenerators := null;
    
    while numGenerators===null and currTrial<maxTrials do
    (
        try
        (
            bblog.debug(" deduceNumGenerators: enter try block ") ;
            
            tmppoint := matrix random(rng^1,rng^(blackBox.numVariables) );
            
            bblog.debug(" deduceNumGenerators:computed tmppoint ") ;
            
            valuesMatrix := blackBox.valuesAt( tmppoint );
            
            bblog.debug(" valuesMatrix computed ") ;
            
            --print valuesMatrix;
            
            assert (numRows valuesMatrix==1);
            numGenerators = numColumns valuesMatrix;            
        )
        then ( computed=true ) 
        else ();
        
        
        currTrial = currTrial+1;
    );
    if not computed then error " failed to deduce number of generators/equations";
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


-- dropDegreeInfo():
-- 
--   for some matrix operations (which ones?), degree information needs to be dropped,
--   which is done by this method. Used in '.jacobianAt' which in turn is  used in the 'padicLift' package.
--
dropDegreeInfo = method();

dropDegreeInfo (Matrix) := Matrix=> (mat)->
(
    return map( (ring mat)^(numRows mat) ,(ring mat)^(numColumns mat), mat );
);


-- deduceJacobianAt()
--
-- Constructs a jacobian at a point supposing that the blackBox implements evaluation at a point.
-- Currently expects that blackBox implements 'valuesAt' and knows number of generators/equations(numGenerators).
-- The second dependency could be removed, as the column number  of the returned 
-- 'valuesAt' evaluation ( a row vector) should be the same as number of generators/equations.
-- 
-- Remark: this stuff is very hard to debug because of the try clauses in the black box.
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

    -- attention, duplicate code!!!
    if (blackBox.withChecks) then
    (
        if (not ( blackBox.valuesAt( point )==0))  then  
        (
            --error("point does not belong to the ideal ! ");
            throw new PointNotOnBlackBox from {"errorMessage" => "deduceJacobianAt: point " | toString point | "does not belong to the object! "}
        )
    );
    -- attention, duplicate code!!!
    
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
            coordinateValue := last coefficients (sub(valueVec_(0,equationIdx), epsRng ), Monomials=>{1  , eps } );
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
    
    blackBoxDummy.withChecks = true;

    computedJac := deduceJacobianAt( blackBoxDummy, point);

    targetJac := dropDegreeInfo sub(jacobian generators I, point); 
    assert( computedJac == targetJac );

    point = matrix{{ 0_R, 0_R, 5_R}};
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

BlackBoxParameterSpace = new Type of  MutableHashTable;




-- question: what is a (succeeded) jet of length 1 at a singular point??






-- jetAtWithInfoResultFunction
-- 
-- constructs a jetAtWithInfoResult  Hashtable from (bestJet and failedJetLength) 
-- with entries "jet", "bestJet", "failedJetLength", "succeeded" 
--
-- purpose: return multiple variables at once in jetAtWithInfo. 
-- Returning a sequence is not an option since then the implementation is not extensible
-- (how to add an additional returned variable without breaking existing code?)
--
-- "bestJet" is always set to bestJet and "failedJetLength" always to failedJetLength.
-- If failedJetLength is null, then "jet" is also set to bestJet and "succeed" to true.
-- otherwise "jet" is set to null and "succeed" to false. 
--
jetAtWithInfoResultFunction := (bestJet, failedJetLength)->
(    
    assert(class bestJet === Jet);
    jet := null;
    if (failedJetLength === null) then 
    (
        jet = bestJet;
    )
    else
    (
        assert(class failedJetLength===ZZ);
        assert(  failedJetLength>=0);
    );
    
    return new HashTable from {    "jet"            => jet, 
                                   "bestJet"        => bestJet, 
                                   "failedJetLength"=> failedJetLength,
                                   "succeeded"      => failedJetLength === null
                              };
);

--
-- see jetAtWithInfoResultFunction()
--
jetAtWithInfoResult := method();
jetAtWithInfoResult ( Jet, ZZ ) := HashTable => (bestJet, failedJetLength)->
(    
    return jetAtWithInfoResultFunction(bestJet, failedJetLength);
);

-- special case for null...
--
-- see jetAtWithInfoResultFunction()
--
jetAtWithInfoResult ( Jet, Nothing ) := HashTable => (bestJet, failedJetLength)->
(
     return jetAtWithInfoResultFunction(bestJet, failedJetLength);
);



-- continueJetWithInfo()
--
-- Ccontinues a given jet up to a requested jetLength (if possible)
-- Returns a hashtable with a bunch of information, see jetAtWithInfoResult()
--
-- todo: eventually cache jacobian and jacobian kernel
--
continueJetWithInfo = method();

continueJetWithInfo( BlackBoxParameterSpace, Jet, ZZ ) := HashTable => ( blackBox,  jet, jetLength )  ->
(  
    assert ( 0 == blackBox.valuesAt(jet#"value") );
    assert ( jetLength >= 0 );
    assert ( jetLength >= length jet );  
    
    failedJetLength := null;    

    if (jetLength==0) then return   jetAtWithInfoResult(jet, failedJetLength);
    
    point := jet#"point";
    
    epsPrecision := length jet;  
          
    coeffRng := (blackBox.coefficientRing); -- we need the braces here !!!
    
    jetValue := jet#"value";
    
    liftingFailed := false;
       
    jetObj := null;
    prejet := null;

    succeededJetLength := length jet;    
    jacobianM2Transposed := transpose blackBox.jacobianAt(point) ;
    
    
    jacobianKernel := generators kernel jacobianM2Transposed ; -- syz would also work
    
    if (length jet==0) then 
    (
        epsPrecision = 1;
        epsRng := getEpsRing( blackBox.coefficientRing,  epsPrecision );
        eps := (gens epsRng)#0;

           
        rnd := random( coeffRng^(numColumns(jacobianKernel)), coeffRng^epsPrecision );
        if (numColumns(jacobianKernel)>0) then 
        (   
            while  zero(rnd) do
            (
                rnd = random( coeffRng^(numColumns(jacobianKernel)), coeffRng^epsPrecision );
            );
        );

        lengthOneLift := sub(point,epsRng) + transpose(sub( jacobianKernel*rnd, epsRng) *eps);
        
        -- first lift will always succeed!
        if ( blackBox.valuesAt(lengthOneLift)!=0 ) then 
        (      
            liftingFailed = true;
            failedJetLength = epsPrecision;
        )
        else
        (
            succeededJetLength = 1;
            jetValue = lengthOneLift;
        );
    );
 
    if (not liftingFailed) then 
    (
        for  epsPrecision in (1 + succeededJetLength)..jetLength do 
        (
            epsRng := getEpsRing( coeffRng, epsPrecision);
            eps := (gens epsRng)#0;
    
            prejet =  sub(jetValue,epsRng);

            valuesAtJet := blackBox.valuesAt(prejet );

            rightHandSide := matrix mutableMatrix( coeffRng, numColumns valuesAtJet ,1 );
            
            if not zero(valuesAtJet) then 
            (           
                rightHandSide = transpose last coefficients (valuesAtJet, Monomials=>{ eps^epsPrecision });
                -- one could also use contract since eps^epsPrec is the highest possible degree
            );
    
            rightHandSide = sub(rightHandSide,coeffRng);
        
            if not (0==rightHandSide % jacobianM2Transposed ) then 
            (
                failedJetLength = epsPrecision;
                liftingFailed = true;
                break; 
            );
            succeededJetLength = epsPrecision;
            x := rightHandSide // jacobianM2Transposed ;
            x = x + jacobianKernel* random(coeffRng^(numColumns(jacobianKernel)), coeffRng^1 );
            x = transpose x;
    
            nextJetValue := sub (prejet, epsRng ) - sub( x, epsRng ) * eps^epsPrecision;
            assert ( 0 == blackBox.valuesAt(nextJetValue) ); -- debug
            jetValue = nextJetValue;
        );
    );

    bestJetObject := jetObject (blackBox,  point, jetValue, succeededJetLength);
    
    return  jetAtWithInfoResult(bestJetObject, failedJetLength);
);

TEST ///
  -- test for bug that  valuesAt(jet) is non zero (jet was incorrectly in a ring with higher precision than required)

kk = QQ
R = QQ[x,y]
I = ideal (x*y)

origin = matrix{{0,0_kk}}
p1     = matrix{{1,0_kk}}

myValuesAt = (p) -> (  return gens sub(I,p);  );

bb = blackBoxIdealFromEvaluation(2,kk, myValuesAt)

jetStats = jetStatsAt(bb,origin,2,10)

jetsL1 = jetStats#"jetSets"#1
jetL1 = first jetsL1#"jets"
L1values =  bb.valuesAt(jetL1#"value")
assert (0 == L1values)
epsRng = ring L1values;
assert (epsRng#"epsPrecision"==1)
-- check that we are in R[eps]/eps^2 for precision 1 jet
assert (ideal epsRng == ideal (first gens last epsRng.baseRings)^2)

jetL2 = jetAt(bb,p1,2)
valuesL2 =  bb.valuesAt(jetL2#"value")
assert (0 == valuesL2)
epsRng = ring valuesL2;
assert (epsRng#"epsPrecision"==2)
-- check that we are in R[eps]/eps^3 for precision 2 jet
assert (ideal epsRng == ideal (first gens last epsRng.baseRings)^3)

jetL0 = jetAt(bb,p1,0)
assert (jetL0#"value" == p1)

///

-- jetAtWithInfo(): 
--
--   tries once to compute a jet , ( see http://en.wikipedia.org/wiki/Jet_%28mathematics%29 for jet definition;) 
--   for the used computation algorithm see the bacherlor thesis at 'http://www.centerfocus.de/theses/js.pdf' .
--
--   preconditions: black box provides evaluation at a point ('valuesAt') and valuesAt(point) evaluates to zero.
-- 
--   returns a hashtable with entries
-- 
--  - "succeeded" a boolean, 
--  - "failedJetLength"  contains the jet length at which the computation failed, otherwise null
--  - "jet"  contains the jet, if succeeded, otherwise null. The jet of the length n has the form
--           j = point + eps*y_1 + . . . + eps^n * y_n 
--           such that F(j) = 0, 
--           where F: E_(n+1)^m -> E_(n+1)^k 
--           with E_(n+1) = K[eps]/( eps^(n+1) ) 
--           whereby K is the coefficient ring (blackBox.coefficientRing), 
--           m is the number of variables (blackBox.numVariables) of the parameter space (same as entries in the point vector)
--           and k is the number of the generators/equation of the (implicitly or explicitly) given ideal. 
--

 
jetAtWithInfo = method();

-- jetAtWithInfo():
--
-- here we improve precision by 1 in each step
-- using Newtons-Algorithm one could double precision in each step, but
-- for this we also need high precision jacobi-matrices.
-- For black-Box-Jacobi-Matrices we might not have high precision Jacobi-Matrices
-- todo question: what do we mean by high precision Jacobi-Matrices?
-- todo: remove duplicate code (see continueJetWithInfo)
--
jetAtWithInfo( BlackBoxParameterSpace, Matrix, ZZ ) := HashTable => ( blackBox,  point, jetLength )  ->
(
    assert ( jetLength >= 0 );
    
    if not (blackBox.isZeroAt(point)) then 
    (
        --error(" point is not on BlackBox ");
        throw new PointNotOnBlackBox from {"errorMessage" => "jetAtWithInfo: point " | toString point | "does not belong to the object! "}
    );

    liftingFailed := false;
    failedJetLength := null;
    jet := point;

    succeededJetLength := 0;    
    
    
    localJetObject := jetObject (blackBox,  point, jet, succeededJetLength);
    
    if (jetLength==0) then 
    (
        return  jetAtWithInfoResult(localJetObject, failedJetLength);
    );
    
    return continueJetWithInfo(blackBox, localJetObject, jetLength);
);





-- jetAtOrException ()
--
-- Computes a jet with given jetLength once using jetAtWithInfo()
-- Returns the jet if succeeded, otherwise the point is singular and an SingularPointException is thrown
-- 
jetAtOrException = method();

jetAtOrException( BlackBoxParameterSpace, Matrix, ZZ) := Jet => ( blackBox,  point, jetLength )  ->
(
    jetResult  := jetAtWithInfo ( blackBox,  point, jetLength);      
    
    if (jetResult#"jet"=== null) then 
    (
       --error ("point is not smooth"); -- is better 
       throw new SingularPointException from {"errorMessage"=>"Point is not smooth",
                                              "failedJetLength" => jetResult#"failedJetLength",
                                              "failedJet" => jetResult#"bestJet",                                            
                                                };
    );     
    return jetResult#"jet";
);



-- continueJetOrException()
--
-- Continues a given jet up to a requested jetLength (if possible) using continueJetWithInfo()
-- returns the computed jet if succeeded, otherwise the point is singular and a SingularPointException is thrown.
--
--
continueJetOrException = method();

continueJetOrException( BlackBoxParameterSpace, Jet, ZZ) := Jet => ( blackBox,  jet, jetLength )  ->
(
    jetResult  := continueJetWithInfo ( blackBox,  jet, jetLength);      
    
    if (jetResult#"jet"=== null) then 
    (
       throw new SingularPointException from {"errorMessage"=>"Point is not smooth",
                                              "failedJetLength" => jetResult#"failedJetLength",
                                              "failedJet" => jetResult#"bestJet",                                            
                                                };
    );     
    return jetResult#"jet";
);
 


-- jetAtOrNull ()
--
-- computes a jet with given jetLength once using jetAtWithInfo()
-- returns the computed jet if succeeded, otherwise returns null.
-- 

jetAtOrNull = method();
jetAtOrNull( BlackBoxParameterSpace, Matrix, ZZ) := Jet => ( blackBox,  point, jetLength )  ->
(
    
    jetResult  := jetAtWithInfo ( blackBox,  point, jetLength);      
    
    return jetResult#"jet";
);


-- continueJetOrException()
--
-- Continues a given jet up to a requested jetLength (if possible) using continueJetWithInfo()
-- returns the computed jet if succeeded, otherwise null.
--
--
continueJetOrNull = method();
continueJetOrNull( BlackBoxParameterSpace, Jet, ZZ) := Jet => ( blackBox,  jet, jetLength )  ->
(
    
    jetResult  := continueJetWithInfo ( blackBox,  jet, jetLength);      
    
    return jetResult#"jet";
);


-- jetAt()
--
-- is the default for computing a jet at a point. 
--
-- the default is set to jetAtOrException(), because computation at a singular point should hurt.
--
jetAt = jetAtOrException;  


-- continueJet()
--
-- is the default method for continuing a jet at a point. 
--
-- the default is set to continueJetOrException(), because computation at a singular point should hurt.
--
continueJet = continueJetOrException;


-- jetStatsAt()
--
-- computes jet statistics at a point, namely the counts of first failed jet lenght for several trials
-- This may be of interest at singular points.
--
-- Parameters:  a black box, a point, the maximal jet length and number of trials to compute a jet.
--
-- Returns a hashtable with 
--
-- "targetJetLength"
-- "numTrials" 
-- "jetSets" -- a hashtable of jet list at point with their length as HashTable key
-- "failedLengthCount" -- a (Tally) where the key is the jetLength l,
--                      and the value is the count of trials where the computation failed at length l.
--
jetStatsAt = method();
jetStatsAt( BlackBoxParameterSpace, Matrix, ZZ, ZZ) := HashTable => ( blackBox,  point, jetLength, numTrials )  ->
(
   if ( numTrials<1 ) then error "jetAtStats: expected numTrials >=1 ";
    
    jetStats := new MutableHashTable ;
        
    jetStats#"jetSets" = new MutableHashTable ;       
    jetStats#"targetJetLength" = jetLength;
    jetStats#"numTrials" = numTrials;            
    jetStats#"failedLengthCount" = new Tally;
    
    jetStats#"succeededCount" = 0;
        
    jetResult  := null;           
    
    for i in 1..numTrials do
    (
        jetResult = jetAtWithInfo ( blackBox,  point, jetLength);
        if (jetResult#"jet"=== null) then   
        (
            jetStats#"failedLengthCount" = jetStats#"failedLengthCount" + tally {jetResult#"failedJetLength"};
        )
        else
        (
            jetStats#"succeededCount" = jetStats#"succeededCount"+1;
                
        );
        if (jetResult#"bestJet"=!= null) then   
        (
            currentJet := jetResult#"bestJet";
            currentJetLength  := length currentJet;           
            
            if (not (jetStats#"jetSets")#?currentJetLength) then 
            (
                (jetStats#"jetSets")#currentJetLength = new JetSet from currentJet;
            )
            else
            (
                addElement(jetStats#"jetSets"#currentJetLength, currentJet);
            );                                    
        );
    );
    jetStats#"jetSets" = new HashTable from    jetStats#"jetSets";
    return new HashTable from jetStats;
);


--
--
--
isCertainlySingularAt = method();
isCertainlySingularAt( BlackBoxParameterSpace, Matrix, ZZ, ZZ) := MutableHashTable => ( blackBox,  point, jetLength, numTrials ) ->
(
    if ( numTrials<1 ) then error "isCertainlySingularAt: expected numTrials >=1 ";
    
    for i in 1..numTrials do
    (
        jetOrError := catch jetAt( blackBox,  point, jetLength);
        if isDerivedFrom(jetOrError,SingularPointException) then 
        (
            return true;
        );
        if not isDerivedFrom(jetOrError, Jet) then 
        (
            throw jetOrError; -- it is a different error 
        );
    );
    return false;
);




isCertainlySingularAt( BlackBoxParameterSpace, Matrix, HashTable) := MutableHashTable => ( blackBox,  point, options ) ->
(
    return isCertainlySingularAt(blackBox, point, options.precision, options.numTrials);
);


-- keysWithoutSymbols(): 
-- 
-- returns hashtable keys without symbol keys.
--
keysWithoutSymbols = method();
keysWithoutSymbols(HashTable) := List => (bb)->
(
    bbh:= new HashTable from bb;
    keylist := keys bbh;
    keylistResult := keylist;
    for key in keylist do 
    (  
        if (class key)===Symbol then
        (
            keyAsString := toString key;
            if bbh#?keyAsString then
                keylistResult=delete(key,keylistResult);
        );
    );
    return keylistResult;
)



--
-- disable BlackBoxParameterSpace construction using new for an arbitrary parameter:
--
new BlackBoxParameterSpace from Thing := ( E, thing) -> 
(
    error "creating blackbox from  type " | toString E | " not implemented ";
);

--
-- disable BlackBoxParameterSpace construction using new without parameters:
--
new BlackBoxParameterSpace   := (E) -> 
(
    error "creating empty blackbox not possible. You have at least to provide the number of variables and their ring";
);


BlackBoxIdeal = new Type of  BlackBoxParameterSpace;


-- disable BlackBoxIdeal construction using new without parameters ( bb = new BlackBoxIdeal ):
new BlackBoxIdeal   := (E) -> 
(
    error "creating empty blackbox not possible. You have at least to provide the number of variables and their ring";
);

-- disable BlackBoxIdeal construction using new for an arbitrary parameter. ( bb=new BlackBoxIdeal from ... ):
new BlackBoxIdeal from Thing := ( E, thing) -> 
(
    error "creating blackbox from  type " | toString E | " not implemented ";
);




-- Type MapHelper contains a matrix and a function
-- to evaluate this matrix on a point
-- This is for example used when projecting onto a
-- subspace (i.e elimination of variables)
--
MapHelper = new Type of HashTable;

--
--
--
createMapHelper = (mapMatrix, imageRing) -> 
(
    mapData := new MutableHashTable;
    
    mapData#"imageRing" = imageRing;
    mapData#"matrix" = mapMatrix;
    
    mapData#"valueAt" =  method();    
    mapData#"valueAt" (Matrix) := Matrix => (point)->
    (
        return sub(mapMatrix,point);
    );
   
    mapData#"valueAtJet" = method();
    mapData#"valueAtJet" (Jet) := Jet => (jet) -> 
    (
        return jetObject(jet#"parent", jet#"point", (mapData#"valueAt")(jet#"value"), jet#"jetLength");                  
    );
   
    return new MapHelper from mapData
);

new MapHelper from Matrix := (XXX, mapMatrix) -> 
(
    -- das hier ist irgendwie alles Quatsch...
    --sourceRing := ring mapMatrix;
    --K := coefficientRing sourceRing;
    --m := rank source mapMatrix;
    --xxx := symbol xxx;    -- todo: get symbol in user space?
    --imageRing := K[xxx_1..xxx_m];
    imageRing := null; 
    return createMapHelper(mapMatrix, imageRing);
);



JetSet = new Type of MutableHashTable;
 

new JetSet from Thing := ( E, thing) -> 
(
    error "creating JetSet from  type " | toString E | " not implemented or not possible ";
);


new JetSet from Jet := ( E, thing) -> 
(
   mht := new MutableHashTable;
   
   mht#"jets" = new List from {thing};
   mht#"point" = thing#"point";
   return mht;
);

net (JetSet) := Net => (jetSet)->
(
    result := net "JetSet{ point: " |net jetSet#"point" | net ", jets{.."| net size jetSet | net "..}}";
    return result;
);

firstElem = method();

firstElem (JetSet) := Jet => (jetset)->
(
    return first jetset#"jets";
);

maximalConditions = method();
maximalConditions( JetSet ) := ZZ => (jetSet)->
(
    return 1+sum apply(jetSet#"jets", jet->length jet);
)


-- how to hide compatible?
compatible = method();


compatible (Jet,Jet) := Boolean => (jet1, jet2 )->
(    
    if  ( jet1#"parent"===jet2#"parent" and 
          jet1#"point" ===jet2#"point"      ) then
         (
            return true;
         );
         return false;        
);

compatible (JetSet,JetSet) := Boolean => (jetSet1, jetSet2 )->
(    
    if (0== size jetSet1 or 0== size jetSet2 ) then   return true;   
    return compatible(firstElem jetSet1, firstElem jetSet2);
);

size (JetSet) := ZZ =>(jetset)->
(
    return #(jetset#"jets");
)


-- probablySameComponent; certainlyDifferentComponent

addElement = method();

addElement (JetSet,Jet) := JetSet => (jetSet, jet)->
(
    if ( size jetSet===0 or
         compatible(firstElem jetSet, jet) 
       ) then
        (
            jetSet#"jets" = append(jetSet#"jets",jet);            
            return jetSet;
        );
    error ("JetSet and Jet are probably incompatible");
)

joinJetSets = method();

joinJetSets (JetSet,JetSet) := JetSet => (jetSet1, jetSet2 )->
(
    if (compatible (jetSet1, jetSet2)) then 
    (
        result := new JetSet;
        result#"jets" = unique join (jetSet1#"jets",jetSet2#"jets");
        return result;
    )
    else
    (
        error ("jet sets probably not compatible");
    );
);




TEST ///
    -- restart 
    -- loadPackage "BlackBoxIdeals"
    debug BlackBoxIdeals --otherwise we have no 'compatible()'
    K = QQ;
    R = K[x,y];
    bbI = new BlackBoxIdeal from ideal x*y;
    
    
    point = matrix{{0_QQ, 0_QQ}};
    
    jet =  jetAt(bbI,point,1)
    
    --locus jet
    jet2 =  jetAt(bbI,point,1)   
    
    compatible(jet,jet2) 
    
    jetSet = new JetSet from jet
    
    --locus jetSet      
    compatible(jetSet,jetSet)    
    
    addElement(jetSet,jet2)
    
    joinJetSets(jetSet,jetSet);    
    -- unique: not im    
///


load "./BlackBoxIdeals/Interpolation.m2";


-- development remark: we need jetAt at least per blackbox individually.

pointProperties = method();

pointProperties (BlackBoxParameterSpace) := List => (bb)->
{
    return bb.pointProperties();
}


attributes = method();

attributes (BlackBoxParameterSpace) := List => (bb)->
{
    return bb.attributes();
}


memberMethods = method();

memberMethods (BlackBoxParameterSpace) := List => (bb)->
{
    return bb.memberMethods();
}

-- todo: we have to differ between interpolator and interpolation.


-- internal method to create a basic black box ( a black box for a parameter space )
--
-- The reason for using an internal (mutable black box) and a public protected blackbox 
-- is to prevent the user from accidental object changing.

-- The object write protection has to be used as  follows: 
--  as a final object for the user always a nonmutable object copy has to be returned;
--  while for the internal 'class' inheritance the non-copied original mutable object 
--  (as created by the internal methods) is needed.
--
-- Since the user has access only to the copy and may modify the object e.g. by registering properties,
-- all (potential) mutable black box properties needs to be stored as local variables in the internal methods, 
-- like 'bbPointProperties' in 'blackBoxParameterSpaceInternal'.
--
-- Access to potential mutable  black box properties by a user can only be modelled  through 'get'-methods, 
-- like blackBox.numGenerators()
-- and may never accessed directly, because then other (shallow) copies of the same black box would run out of sync.
-- Because blackBox.numGenerators() is defined in the same context with the 
-- local variable 'localNumGenerators', the variable 'localNumGenerators'
-- is visible (and modifiyable) inside 'blackBoxParameterSpaceInternal', but not outside!


blackBoxParameterSpaceInternal = method();

blackBoxParameterSpaceInternal( Type, ZZ, Ring  ) := HashTable => ( resultType, numVariables, coeffRing ) ->
(
    
    blackBox := new MutableHashTable;
    
    blackBox = newClass(resultType, blackBox);
    
    blackBox.withChecks = true;
    
    blackBox.disableChecks = ()->
    (
         blackBox.withChecks = false;
    );
    
    blackBox.enableChecks =  ()->
    (
         blackBox.withChecks = true;
    );
    
    -- public: 
    blackBox.coefficientRing = coeffRing;  
    blackBox.numVariables = numVariables;  -- stores the number of the variables in the parameter space.
    
    -- private: 
    bbPointProperties := new HashTable  ;     -- stores the point properties
    localNumGenerators := null;                  -- stores the number of the generators/equations. Since this value may be updated,
                                           -- the value cannot be stored in the blackBox directly as a key-value pair, 
                                           -- because, otherwise different black box references referring to the same object
                                           -- could get out of sync. The variable is accessed by a getter(numGenerators())
    
    
    
    
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
            try ( if coefficientRing ring point=!= blackBox.coefficientRing then error errorMsg )  
            then () 
            else ( error (errorMsg) );
        );
        if (numColumns point=!=blackBox.numVariables) then 
        ( 
            error (" illegal point : expected " | toString blackBox.numVariables | " coordinates");
        );
    );


    -- pointProperty():
    --
    --    return a point property (by symbol or by name) stored in 'bbPointProperties'. 
    --    If a property does not exist, throws an error.
    --
    blackBox.pointProperty = method();

    --
    -- get a point property by name
    --
    blackBox.pointProperty (String) := Function =>( prop ) ->
    (
        if not bbPointProperties#?prop then        
            error (" blackbox does not have requested point property " | toString prop  );

        return bbPointProperties#prop;
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
        return bbPointProperties#?propertyName;
    );
    
    blackBox.hasPointProperty ( Symbol ) := Boolean => (propertySymbol)->
    (
        return blackBox.hasPointProperty(toString propertySymbol)
    );

    
    -- 'pointProperties':   
    --     returns a list of all registered point properties (as Strings)
    --
    blackBox.pointProperties = ()->
    (   
        return sort unique apply (keys bbPointProperties, key->toString key);
    );


    -- 'pointPropertiesAsSymbols' : 
    --  returns a list of all registered point properties (as Symbols)
    --

    getPropertySymbol := method ();
    getPropertySymbol(String) := Symbol=> (propertyName)->
    (
     
        propertySymbol := null;
        try  (  propertySymbol = getGlobalSymbol propertyName; ) 
        else 
        ( 
            propertySymbol = getGlobalSymbol( User#"private dictionary", propertyName); 
        );
        return propertySymbol;
    );

    blackBox.pointPropertiesAsSymbols = ()->
    (   
        return apply( blackBox.pointProperties(), propertyName-> getPropertySymbol(propertyName) );
    );


    -- setPointProperty(): 
    --  
    --  internal method to set a point property. 
    --                     Is called by 'outerSetPointProperty', 'setIsZeroAt', 'setValuesAt', 'setJacobianAt'
    -- 
    -- the method works as follows: 
    -- 1. the current (internal) variable 'bbPointProperties' is copied and transformed to a mutable HashTable
    -- 3. the propertyMethod (see 2) is added to the internal variable  'bbPointProperties'
    -- 4. the 'bbPointProperties' are changed to immutable.
    -- 5. the (internal) blackBox HashTable is extended by methods which accept one parameter and
    --     call the corresponding propertyMethod in 'pointProperty' This level of indirection is done,
    --     to keep access to the correct bbPointProperties even if the internal BlackBox object is replaced or updated: 
    --     all (updated) blackboxes will refer to the same 'bbPointProperties' variable. 
    --  

    --  after a call of 'getUpdatedBlackBox', the property is accessible by 
    --  its symbol name and eventually by the symbol, if there is no symbol clash
    -- 
    -- Remark. the first parameter is a symbol and not a name, because it is imaginable, that a user / a package author 
    --   could want to pass the concrete symbol he wants.  
    --
    setPointProperty := method();
    setPointProperty ( Symbol, Function ) := Thing => ( propertySymbol, propertyMethod )->
    (

        bblog.debug(" called setPointProperty for " | toString propertySymbol ) ;

        propertyName := toString propertySymbol;

        assert(propertyName=!=null); 

        packageSymbol := null;
        if BlackBoxIdeals.Dictionary#?propertyName then 
        packageSymbol  = BlackBoxIdeals.Dictionary#propertyName;

        -- step 1
        bbPointProperties = new MutableHashTable from bbPointProperties;

        bblog.debug(" set point Property " | propertyName) ;

 
        -- step 2,3
        -- remember: cannot return in a try statement!
        bbPointProperties#propertyName = ( point )-> 
        ( 
            checkInputPoint(point);
            result := (propertyMethod)(  point ); 
            return result;
        );

        -- inconsistent and unnecessary ?    
        -- bbPointProperties#propertySymbol = bbPointProperties#propertyName;
    
        if packageSymbol=!=null then 
        (
            bbPointProperties#packageSymbol = bbPointProperties#propertyName;
        );
  
        -- step 5
        if mutable blackBox then 
        (
            bblog.debug(" mutable blackBox ") ;
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
        bbPointProperties = new HashTable from bbPointProperties;
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
        try  ( assert(numRows result==1); ) else { error ( "valuesAt did not return an one-row matrix "); };
        return result;
    );

    -- setIsZeroAt
    -- 
    --   sets the check, if a point belongs to the object (isZeroAt(point)==true ) or not.
    --
    --   called by   'outerSetPointProperty'<-{'registerPointProperty', 'updatePointProperty'}, 'setValuesAt'
    -- 
    setIsZeroAt := (pIsZeroAt) ->
    (   
        -- parameter is called differently to the symbol 'isZeroAt', otherwise it seems we could get the wrong value...
        setPointProperty("isZeroAt" , pIsZeroAt );
    );
   
    localJetAt := jetAt;

    
    -- setJacobianAt():
    -- 
    --    set a method to compute the jacobian at a point.
    --
    --   called by 'outerSetPointProperty'<-{'registerPointProperty', 'updatePointProperty'}, 'setValuesAt'
    --  
    --   triggers updates for  isProbablySmoothAt, 'rankJacobianAt'.
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
            return isCertainlySingularAt( blackBox,  point, singularTestOptions );
        );

        setPointProperty("isCertainlySingularAt" , localIsCertainlySingularAtWrapper );

        localIsProbablySmoothAtWrapper  := ( point )  ->
        (
            return not isCertainlySingularAt( blackBox,  point, singularTestOptions );
        );

        setPointProperty("isProbablySmoothAt" , localIsProbablySmoothAtWrapper );


        setPointProperty( "rankJacobianAt" , localRankJacobianAt );
    
        blackBox#"jetAt" = (point, jetLength)->
        (
            return localJetAt(blackBox, point, jetLength);
        );
        
        blackBox.jetAt = blackBox#"jetAt";        
    );

    
    
    updateSingularityTest := ()->
    (
        if blackBox.hasPointProperty("valuesAt") then
        (
            localIsCertainlySingularWrapper  := (  point )  ->
            (
                return isCertainlySingularAt( blackBox,  point, singularTestOptions );
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
                return not isCertainlySingularAt( blackBox,  point, singularTestOptions );
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
        return;
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

        blackBox.numGenerators = ()->
        (
            -- if (localNumGenerators===null) then ( error "failed to deduce number of generating equations" );
            return localNumGenerators;
        );

        bblog.info( "updated blackBox.numGenerators to " | toString blackBox.numGenerators() );   

        setIsZeroAt(
            (point)->( return blackBox.valuesAt(point)==0 ;) 
        );

        ----- jacobian at:
        localMethod :=  (point)->deduceJacobianAt( blackBox, point );
        setJacobianAt ( localMethod );  
        updateSingularityTest();
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
        --    error (" provided method " | propertyName | " expected to accept 2 parameters ( blackbox, point ),  
        -- but the passed one seems to accept " | toString acceptedNumParameters);
        
        -- now wrap the provided method if neccesary in a way that it accepts only a point: 

        localPropertyMethod := propertyMethod;

        if acceptedNumParameters==2 then 
        (
            localPropertyMethod = ( point )-> 
            ( 
                return propertyMethod( blackBox,   point ); 
            );
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

    blackBox.setPointProperty = method();

    blackBox.setPointProperty(String, Function) := (BlackBoxParameterSpace) => ( propertyName, propertyMethod )->
    (
        propertySymbol := getPropertySymbol(propertyName);
        assert(propertySymbol=!=null);
        outerSetPointProperty( propertySymbol, propertyMethod );
        blackBox.updateBlackBox();
        return blackBox;
    );

    --  updatePointProperty()
    --
    -- if a property is already set, updates it.
    --
    -- todo: test if updating 'valuesAt' will  trigger updating isZeroAt and jacobianAt 
    --     
    blackBox.updatePointProperty = method();

    blackBox.updatePointProperty(String, Function) := (BlackBoxParameterSpace) => ( propertyName, propertyMethod )->
    (
        if  (  bbPointProperties#?propertyName ) then
        (  
            return  blackBox.setPointProperty( propertyName, propertyMethod );
        ) 
        else 
        (
            error ("point property "| toString propertyName | " does not exist.") 
        );
    );
   

    blackBox.updatePointProperty(Function) := (BlackBoxParameterSpace) => (  propertyMethod )->
    (
        propertyName := toString propertyMethod;
        return updatePointProperty(propertyName, propertyMethod);
    );
    
    -- updateBlackBox()
    --
    -- update keys of the blackBox Hashtable in case there are known point properties but no corresponding 
    -- keys. Initial purpose (that changed): blackbox variable was by intention not writeable and modification needed
    -- copying. 
    -- Current purpose: access point property by property name string via '#' operator. 
    -- But, question, does this also add the symbolic stuff??? Something seems still weird here...(jk, 07.02.2017)
    -- 
    blackBox.updateBlackBox = () ->
    (
        return;
        -- not necessary anymore(?)
        for  property in blackBox.pointProperties() do
        (
            propkeys := unique sort {( property )} | {toString property};
           
            for key in propkeys  do
            (
                if not blackBox#?key then 
                (
                    -- hier springen wir jetzt nie(?) rein.
                    --print("here3");
                    blackBox.pointProperty(property);
                    --blackBox#key = (point)->( (blackBox.pointProperty(toString property))(point) );
                    blackBox#key =  (blackBox.pointProperty(toString property)) ;                 
                )
                else 
                (      
                );     
            );    

        );
    );


    -- registerPointProperty()
    --
    -- a method to register a point property, while providing a propertySymbol,  
    -- expecting that after registering (and getUpdatedBlackBox() ) the property will be accessible via  blackBox#propertySymbol .
    -- Usually providing the corresponding symbol is not necessary, but it could be, since each package has its own symbol scope.
    --
  
    --
    
    blackBox.registerPointProperty = method();    
    
    -- todo: this method with this interface should not be publicly visible or accessible (is an internal one)
    --
    blackBox.registerPointProperty(String, Symbol, Function) := BlackBoxParameterSpace => 
      ( propertyName, propertySymbol, propertyMethod )->
    (
        assert( (toString propertySymbol)==propertyName);

        if  ( not  bbPointProperties#?propertyName 
        and not  blackBox#?propertySymbol          and   not  blackBox#?propertyName  ) then 
        (
            outerSetPointProperty( propertySymbol, propertyMethod );
        )
        else error(" key  "| propertyName |"  exists already. If it is a point property, please use 'updatePointProperty' for updating.");
        blackBox.updateBlackBox();
        return blackBox;
    );

    -- public..
    --
    blackBox.registerPointProperty(String, Function) := Thing => ( propertyName, propertyMethod )->
    (
        propertySymbol :=  getPropertySymbol(propertyName);
        return blackBox.registerPointProperty(  propertyName, propertySymbol, propertyMethod )
    );
    
    blackBox.registerPointProperty( Function) := Thing => (  propertyMethod )->
    (
        propertyName := toString propertyMethod;
        propertySymbol :=  getPropertySymbol(propertyName);
        return blackBox.registerPointProperty(  propertyName, propertySymbol, propertyMethod )
    );

    blackBox.rpp =  blackBox.registerPointProperty;

    blackBox.upp =  blackBox.updatePointProperty;

    --
    -- memberMethods()
    --
    -- return a list of known methods. Manually updated.
    --
    blackBox.memberMethods = ()->
    (   
        methods:= { getGlobalSymbol( BlackBoxIdeals.Dictionary, "memberMethods" ) ,
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "attributes" ) ,
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "pointProperties" ),
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "pointPropertiesAsSymbols" ),
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "hasPointProperty" ),
                    --getGlobalSymbol( BlackBoxIdeals.Dictionary, "pointProperty" ),
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "registerPointProperty" ),
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "setSingularityTestOptions" ),
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "rpp" ), 
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "upp" ), 
                    getGlobalSymbol( BlackBoxIdeals.Dictionary, "updatePointProperty" )
                    --getGlobalSymbol( BlackBoxIdeals.Dictionary, "unknownIsValid" )
            };
    --  methods:= {   memberMethods,
    --                attributes,
    --                pointProperties,
    --                pointPropertiesAsSymbols,
    --                hasPointProperty,
    --                pointProperty,
    --                registerPointProperty, 
    --                updatePointProperty,
    --                unknownIsValid
    --              };

        sortedMethods := sort apply(methods, i-> ( toString i, i ));
        sortedMethods = apply(sortedMethods, i->(i_1));
        return sortedMethods;  
        --return methods;
    );
  
  
    -- attributes()
    --
    -- returns a list of known attributes. 
    -- (computed as 'keys blackBox which are not Functions' \ {pointPropertiesAsSymbols() }
    --
    blackBox.attributes = ()->
    (
        all :=  keysWithoutSymbols blackBox;
        all = select (all, (foo)->(not instance(blackBox#foo,Function)));   
        kM := kPP := kPS := {};
        -- kM =  blackBox.memberMethods();
        -- kPP =  blackBox.pointProperties();
        kPS  =  blackBox.pointPropertiesAsSymbols();
        toRemove := kM | kPP |kPS;
        for symb in toRemove do
        all = delete(symb,all);
        if ( blackBox#?"valuesAt" ) then
        (
            all = all | { getGlobalSymbol( BlackBoxIdeals.Dictionary, "numGenerators" ) };
        );
        sortedAttributes := sort apply(all, i-> ( toString i, i ));
        sortedAttributes = apply(sortedAttributes, i->(i_1));
        return sortedAttributes;  
    );


    -- singularityTestOptions()
    --
    -- returns currently used configuration for singularity test (at a point)
    --
    blackBox.singularityTestOptions = ()->
    (
        if not blackBox.hasPointProperty("valuesAt") then
        ( 
                error "no singularity test options: valuesAt-property not available";
        );
        return new HashTable from singularTestOptions;
    );


    -- setSingularityTestOptions()
    --
    -- sets currently used configuration for singularity test (at a point)
    --
    blackBox.setSingularityTestOptions = (prec, numTrials)->
    (
        if not blackBox.hasPointProperty("valuesAt") then
        ( 
                error "cannot set singularity test options: valuesAt-property not available";
        );

        if (prec <0) then error "setSingularityTestOptions: expected prec >= 0";
        if (numTrials <=0) then error "setSingularityTestOptions: expected numTrials > 0";

        singularTestOptions.precision = prec;
        singularTestOptions.numTrials = numTrials;
        updateSingularityTest();

    );

    
    -- todo: choose later a different component calculator depending on strategy option
    
    -- better: check component calculator type /interface.
    blackBox.setInterpolator = (interpolatorParam)->
    (
        blackBox.onComponentAnswerStrategies = interpolatorParam.onComponentAnswerStrategies;
        blackBox.setOnComponentAnswerStrategy   =  interpolatorParam.setOnComponentAnswerStrategy;
        blackBox.onComponentAnswerStrategy   =  interpolatorParam.onComponentAnswerStrategy;
        
        -- TODO well, setOnComponentPrecision() should be more generic as different component calculators may 
        -- have different configuration settings ( so may be another component calculator has no 'onComponentPrecision' property)
        blackBox.setOnComponentPrecision   =  interpolatorParam.setOnComponentPrecision;
        blackBox.onComponentPrecision   =  interpolatorParam.onComponentPrecision;
        blackBox.interpolator   = interpolatorParam;
        blackBox.resetInterpolation   = interpolatorParam.resetInterpolation;
        blackBox.isOnInterpolatedComponent  = interpolatorParam.isOnComponent;
        blackBox.interpolatedComponents     =  interpolatorParam.components;
        blackBox.interpolatedComponentNames  =  interpolatorParam.componentNames;
        blackBox.componentNameInUse  =  interpolatorParam.componentNameInUse;
        blackBox.interpolatedComponentByName = interpolatorParam.componentByName;
        blackBox.renameComponent =  interpolatorParam.renameComponent;
         -- todo: when we change the interpolator, this will stop to work:
        blackBox.interpolateComponentsAt  = interpolatorParam.interpolateComponentsAt;    
        blackBox.refineInterpolation  = interpolatorParam.refineInterpolation;    
        blackBox.interpolateComponentAt      =  interpolatorParam.interpolateComponentAt;
        blackBox.interpolatedComponents         =   interpolatorParam.components;  
        
        localinterpolatedComponentsAt := (point) -> interpolatorParam.componentsAt(point);
        
        blackBox.setPointProperty("interpolatedComponentsAt", localinterpolatedComponentsAt );    
        
        localinterpolatedComponentNamesAt := (point) -> interpolatorParam.componentNamesAt(point);
        
        blackBox.setPointProperty("interpolatedComponentNamesAt", localinterpolatedComponentNamesAt);
                
        -- TODO well , setSameComponentPrecision should be more generic as different component calculators may 
        -- have different configuration settings ( so may be another component calculator has no 'onComponentPrecision' property)
        --blackBox.setSameComponentPrecision =  interpolatorParam.setSameComponentPrecision;
    );
    
    blackBox.setInterpolator( createSimpleInterpolator(blackBox) );    
    
    -- a user should not call this method...

    return blackBox;
)

TEST ///
  -- test interpolateComponents
  kk := ZZ
  R := ZZ[x,y]
  I:= ideal (x*y*(x^2-y^2))
  bb = new BlackBoxIdeal from I
  point1 = matrix {{1,1_ZZ}}
  point2 = matrix {{1,0_ZZ}}
  pointList = {point1,point2}
  maxDegree := 4
  bb.interpolateComponentsAt(pointList, maxDegree);
///
 


-- blackBoxParameterSpaceInternal()
-- 
-- this function may be used to create a derived object, which inherits properties of an black box ideal
-- since the blackbox object is not copied 
--
blackBoxParameterSpaceInternal(Type, Ring) := HashTable => ( resultType, pRing ) ->
(
    blackBox := blackBoxParameterSpaceInternal(resultType, #(gens pRing), coefficientRing pRing);
    
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

listToStack := (L)->
(
    return ( apply(L, i -> ("--   " | toString i)));
)


net (BlackBoxParameterSpace) := Net =>(bb)->
(
    L := {"--" | toString class bb};
    L = L | {"--"};
    L = L | {"-- attributes:"};
    L = L | {"-- "} | listToStack bb.attributes();
    L = L | {"--"};
    L = L | {"-- methods:"};
    L = L | {"-- "} | listToStack  bb.memberMethods();
    L = L | {"--"};
    L = L | {"-- point properties:"};
    L = L | {"-- "} | listToStack  bb.pointProperties();

    return stack L;
);


blackBoxIdeal = method();

 
blackBoxParameterSpace = method();

--
-- this function is final, that means nobody should use this method for creating a derived object
--
new BlackBoxParameterSpace from Ring := (E, pRing )->
(
    blackBox := blackBoxParameterSpaceInternal(BlackBoxParameterSpace, pRing);
   
    return blackBox;
)

--
-- this function is final, that means nobody should use this method for creating a derived object
--
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
    blackBox := blackBoxParameterSpaceInternal(BlackBoxParameterSpace, numVariables, coeffRing );
    --bb := newClass( BlackBoxParameterSpace, blackBox );
    return blackBox;
)


-- todo: how to check, if 'ring equationsIdeal' is not a quotient ring?

--
-- this function is final, that means nobody should use this method for creating a derived object
--
new BlackBoxIdeal from Ideal := (E, equationsIdeal)->
(
     blackBox :=  blackBoxParameterSpaceInternal( BlackBoxIdeal, ring equationsIdeal );
    
    
    -- maybe blackBox.addProperty( ideal, equationsIdeal)
    blackBox.ideal = equationsIdeal;      
    
    -- maybe blackBox.addProperty( equations, gens equationsIdeal )
    blackBox.equations =  gens  equationsIdeal; 
    
    
    -- registering valuesAt generates 'isZeroAt' and 'jacobianAt', too !  
    blackBox.registerPointProperty("valuesAt",( bb, point)->  
        (
            result :=  gens sub( equationsIdeal , point);
            return  result;
        )
    );


    -- maybe blackBox.addProperty( jacobian, jacobian gens  equationsIdeal )
    blackBox.jacobian = jacobian gens  equationsIdeal;


    -- we have to call updatePointProperty for "jacobianAt", because "jacobianAt" is present. 

    blackBox.updatePointProperty( "jacobianAt",
        ( bb, point )->
            (  
                -- attention, duplicate code!!!
                if (blackBox.withChecks) then
                (
                    if (not ( blackBox.valuesAt( point )==0))  then  
                    (
                        --error("point does not belong to the ideal ! ");
                        throw new PointNotOnBlackBox from {"errorMessage" => " Point " |toString point|" not on Black box"};
                    );
                );
                -- attention, duplicate code!!!
                jacobianM2MatrixAt := sub( blackBox.jacobian , point);
                return jacobianM2MatrixAt;
            )
    );   
    
    return blackBox; 
)

--
-- this function is final, that means nobody should use this method for creating a derived object
--
blackBoxIdeal (Ideal) := BlackBoxIdeal =>(equationsIdeal)->
(
   return new BlackBoxIdeal from equationsIdeal;
)

blackBoxIdeal (Ideal, Boolean) := BlackBoxIdeal =>(equationsIdeal, withChecks)->
(
    bb := new BlackBoxIdeal from equationsIdeal;
    bb.disableChecks();
    return bb;
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
    )  
    else ();
    point = sub( point, coeffRng ) ;
    try {    jac= IFPBlackBox.jacobianAt(point); } 
    else
    (
        error("testblackBoxIdeal: jacobianAt should succeed  due the coefficient ring of the point matrix matches the ideal coefficient ring ");
    );  
    IFPBlackBox.ring;
    IFPBlackBox.valuesAt(point) ;
    assert(   IFPBlackBox.isZeroAt( point ) );
    assert( IFPBlackBox.jacobianAt(point)==sub( jacobian IFP,point) );
    assert( IFPBlackBox.valuesAt(point)== gens sub(  IFP, point ) );
)



--
--  creates a BlackBoxIdeal from a given evaluation method ('valuesAt') which takes a point (a row matrix)
--
--  parameters: numVariables in the parameter space,
--              coefficientRing , the 
--              valuesAt: a method for evaluating the object at one parameter point.
--
--
-- this function is final, that means nobody should use this method for creating a derived object
--

blackBoxIdealFromEvaluation = method();

blackBoxIdealFromEvaluation(ZZ, Ring, Function) := HashTable => ( numVariables, coeffRing, pValuesAt )  ->
(    
    blackBox := blackBoxParameterSpaceInternal( BlackBoxIdeal, numVariables, coeffRing );

    --sets isZeroAt, jacobianAt, rankJacobianAt and numGenerators
    blackBox.registerPointProperty ("valuesAt", (bb,point)->pValuesAt(point) ); 

    check := ()->
    (
         numVariables :=  blackBox.numVariables;

         point := matrix { apply(numVariables, i-> 0_(blackBox.coefficientRing) ) };
         blackBox.valuesAt( point );
         blackBox.isZeroAt( point );
    );

    check();      
    
    --blackBox = newClass( BlackBoxIdeal, blackBox ); 
    return blackBox ;
)


blackBoxIdealFromEvaluation( Ring, Function ) := HashTable => ( pRing, pValuesAt ) ->
(

   blackBox := blackBoxParameterSpaceInternal(BlackBoxIdeal, pRing );
   blackBox.registerPointProperty ("valuesAt", (bb,point)->pValuesAt(point) ); --sets isZeroAt, jacobianAt, rankJacobianAt and numGenerators

   check := ()->
   (
        numVariables :=  blackBox.numVariables;
        
        point := matrix { apply(numVariables, i-> 0_(blackBox.coefficientRing) ) };
        blackBox.valuesAt( point );
        blackBox.isZeroAt( point );
   );

   check(); 
   --blackBox = newClass( BlackBoxIdeal, blackBox ); 
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
   Headline
        a type
   Description
         Text
            A @TO BlackBoxIdeal @ is a special @TO BlackBoxParameterSpace @.
            It is used when the equations of a stratum of
            a parameter space are known at least implicitly.
            
            \,\, \bullet \, If the equations are known explicitly use 
            @TO blackBoxIdeal @ to create a @TO BlackBoxIdeal @.
            
            \,\, \bullet \, If the equations are known implicitly, i.e.
            an algoithm to evaluate the equations is known,
            then use @TO blackBoxIdealFromEvaluation@
            to create a @TO BlackBoxIdeal @.
            
            In both cases the BlackBoxIdeal has the following
            point properties predefined:

            \,\, \bullet \, @TO "isZeroAt" @ \break
            \,\, \bullet \, @TO "valuesAt" @ \break
            \,\, \bullet \, @TO "jacobianAt" @ \break
            \,\, \bullet \, @TO "rankJacobianAt" @ \break
                        
///

doc ///
   Key
        blackBoxIdealFromEvaluation        
   Headline
        create a  BlackBoxIdeal  from an evaluation method
   Usage   
        blackBoxIdealFromEvaluation( R, evaluationMethod )
   Inputs  
        R:Ring
             that contains the implicitly given Ideal
        evaluationMethod:Function
             accepting a point given as a row matrix 
             (with the same number of entries as the number of
             variables of R)
             and returns 
             the evaluation of the generators of the ideal 
             at this point as a row matrix.
   Outputs
        : BlackBoxIdeal
   Description   
        Text    
           Creates a blackbox describing an implicitly given ideal from an evaluation method \break
           A BlackBoxIdeal is a special @TO BlackBoxParameterSpace @.
           
           This is useful when the generators of an ideal are
           to big to write down, but an algorithm evaluating the generators
           of at a point is known.
           
           A trivial example is the determinant of a n x n matrix. 
           The explicit polynomial has n! terms. For big enough n
           it is impossible to even store the n! terms of it. 
           But for a given point in K^{n*n}, the
           matrix has entries in a field and the Gauss-Algorithm can 
           be used to evaluate the determinant. This takes only O(n^3) steps:
           
           Make an n x n matrix from a vector of length n^2:
        Example
           K = ZZ/11
           n = 10;
           matrixAt = point -> matrix apply(n,i->apply(n,j->point_(i*n+j)_0))              
           testPoint = random(K^1,K^(n*n));          
           matrixAt(testPoint)
        Text
           Now take the determinant. To be able to take this
           as an evaluation method of a BlackBoxIdeal it has to
           take a 1 x n^2 matrix and has to return a 1 x m matrix
           (in this case m=1).
        Example
           detAt = point -> matrix{{det matrixAt(point)}};
           detAt(testPoint)
        Text
           Now make a BlackBox from this evaluation method:
        Example
           R = K[a_(0,0)..a_(n-1,n-1)];
           bbDet = blackBoxIdealFromEvaluation(R,detAt);
        Text
           Evaluate the BlackBox at the test point:
        Example
           bbDet.valuesAt(testPoint) == det(matrixAt(testPoint))
        Text
           Using the Package FiniteFieldExperiments
           heuristic information about the irreducible components
           of a variety defined by an black box
           can be obtained.
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
        blackBoxIdeal(I)
   Inputs  
        I:Ideal
   Outputs
        : BlackBoxIdeal
   Description      
        Text

          Sometimes one wants to use the algorithms implemented
          for black box ideals even though explicit generators
          of the ideal are known. In this case one can create
          a BlackBoxIdeal from an ordinary ideal:
        Example
          K = ZZ/11;
          R = K[x,y,z]
          bb = blackBoxIdeal(ideal(x*y,x*z));
          -- this defines a plane x=0 and a line y=z=0
          bb.isZeroAt(matrix{{0,1,1_K}})          
          bb.isProbablySmoothAt(matrix{{0,1,1_K}})
          bb.isCertainlySingularAt(matrix{{0,0,0_K}})
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
 

    -- assert bbRankM#?(global rankMat); -- fails..
    assert bbRankM#?("rankMat");

    assert( rankMatNew(bbRankM,point) == bbRankM.rankMat(point) );

    assert(bbRankM.coefficientRing===ZZ);

    rankMatNew := (blackBox, point)->4 --also influences bbRankM; because rebuild does not copy; it just exports new registered properties.

    bbRankM.updatePointProperty("rankMat",rankMatNew)

    (bbRankM.pointProperty("rankMat"))(point);

    assert( rankMatNew(bbRankM,point) == bbRankM.rankMat(point) );
    assert( rankMatNew(bbRankM,point) == bbRankM#"rankMat"(point) );
    assert( rankMatNew(bbRankM,point) == (bbRankM.pointProperty("rankMat"))(point) );
    assert( rankMatNew(bbRankM,point) == (bbRankM.pointProperty(getGlobalSymbol "rankMat"))(point) );
    keys bbRankM
    
    valuesAt := ( blackBox, point )-> matrix {{1,2}};

    bbRankM.registerPointProperty( "valuesAt", valuesAt );
    -- that is unfortunate; registering a point property requires a rebuild.   
     
    assert(  bbRankM.hasPointProperty("isZeroAt") );

    assert(  bbRankM.hasPointProperty("jacobianAt") );
 
    assert( bbRankM.numGenerators() =!= null)

    assert( bbRankM.numGenerators() === 2 )
    -- bbRankM.numGenerators()
    
    illegalPoint := matrix {{1,2,3,4,5,6}}; 
   
    try ( bbRankM.rankMat(illegalPoint) ) then ( assert(false) ) else ();

    bbRankM = blackBoxParameterSpace( 5 ,ZZ/7 )

    valuesAt := ( blackBox, point )-> matrix {{1,2}};

    --bbRankM.registerPointProperty("valuesAt",valuesAt);
    bbRankM.rpp("valuesAt",valuesAt);

    point  = sub(point,ZZ/7); 

    bbRankM.valuesAt(point)

    illegalPoint := sub(point,ZZ/2); 

    try ( bbRankM.valuesAt(illegalPoint) ) then ( assert(false) ) else ();

   
///

--load "./BlackBoxIdeals/documentation/blackBoxParameterSpace.m2";


-- how to document 
-- blackBoxParameterSpace(ZZ,Ring) ?



doc ///
    Key
        blackBoxParameterSpace
        BlackBoxParameterSpace
    Headline
        create a BlackBoxParameterSpace.
    Usage   
        blackBoxParameterSpace(d,K)
    Inputs  
        d: ZZ 
         the (affine) dimension of the parameter space
        K: Ring
             a field
    Outputs
        : BlackBoxParameterSpace
    Description
        Text
            A black box parameter space is used to implement parameter spaces
            with their universal families in a pointwise fashion.                      
    
            For a quick start see the @TO "Singularities of cubic surfaces" @-tutorial
            \break \break
            
            A more trivial example is the parameter space of 2x2 matrices
            with entries in a finite field. This is
            an affine 4 dimensional space:
            
        Example
            K = ZZ/11;
            bb = blackBoxParameterSpace(4,K);
        Text
            For each point in this parameter space we have a 2x2-matrix:
        Example
            matrixAt = (point) -> matrix{{point_0_0,point_1_0},{point_2_0,point_3_0}};
            bb = bb.rpp("matrixAt",matrixAt);
            bb.matrixAt(matrix{{1,2,3,4_K}})
        Text
            Properties of the parametrized objects can also be store in
            the Black box. Maybe one is especially interested in the rank
            filtration of the parameter space:
        Example
            rankAt = (point) -> rank matrixAt(point)
            bb = bb.rpp("rankAt",rankAt);
            rankAt(matrix{{1,1,1,1_K}})
            rankAt(matrix{{1,2,3,4_K}})
        Text
            The black box parameter space keeps track of all such properties
            defined:
        Example
            bb.pointProperties()
        Text
            How such properties stratify the parameter space can be systematically
            evaluated by using the package FiniteFieldExperiments.
    Caveat
            There are several special property names: 
            
            If  @TO valuesAt@ is registered or updated, this will
            will implicitly construct   
            @TO isZeroAt@, @TO numGenerators@,  @TO jacobianAt@,
             @TO rankJacobianAt@, @TO isCertainlySingularAt@ and @TO isProbablySmoothAt@.

            This happens automatically if the BlackBox is created using
            @TO blackBoxIdealFromEvaluation@.
            
            If  @TO jacobianAt@ is registered or updated, this will
            will implicitly construct   
             @TO rankJacobianAt@, @TO isCertainlySingularAt@ and @TO isProbablySmoothAt@.
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
            bbC.pointProperties()
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
            bbC.pointProperties()
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
            bbC.pointProperties()
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
            bbC.pointProperties()
        Text
            Calculate the degree of the singular locus for our examples
        Example
            bbC.degreeSingularLocusAt(coeffCubicCone)
            bbC.degreeSingularLocusAt(coeffCubicFermat)
        Text
            Now the BlackBoxParameterSpace hat a number of point properties
        Example
            bbC.pointProperties()
        Text
            These properties can now be used in a finite field experiment
            that studies the statification of our parameter space. Here is a
            simple minded version of such an experiment:
        Example
            tally apply(100,i->bbC.degreeSingularLocusAt(random(K^1,K^20))) 
        Text
            We see that there is an open stratum of smooth cubics. The
            largest closed stratum consists of those cubics with a A1 singularity.
            The package FiniteFieldExperiments helps to do the bookkeeping 
            for such experiments and also provides more detailed interpretation
            of the results.

///

doc ///
    Key
        BlackBoxIdeals
    Headline
          black boxes for implicitly given ideals
    Description
        Text
            The BlackBoxes of this Package come in two flavors:
            
            1) @TO BlackBoxParameterSpace@
            
            This represents a family of algebraic objects over
            an affine space K^n in a pointwise fashion. Together with 
            the package FiniteFieldExperiments this can be
            used to study the stratification of the  parameter
            space K^n with respect to various properties of the algebraic
            objects parametrized. For a quick start look at the
            @TO "Singularities of cubic surfaces" @-tutorial.
            
            2) @TO BlackBoxIdeal@
             
            A BlackBoxIdeal is a special BlackBoxParameterSpace.
            Here equations for a specific stratum of the
            parameter space are known at least implicitly.
            
            Together with 
            the package FiniteFieldExperiments this stratum can
            then be studied much more precisely. For example
            heuristic estimates on the number and codimension
            of its reduced components can be obtained. If one
            is lucky, even equations for the different components
            can be found. For a quick start look at the ?????-Tutorial.

    Caveat
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

    bb = blackBoxIdeal I

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
        singularityTestOptions
    Headline
        show current parameters for singularity test
    Usage   
        bb.singularityTestOptions()
    Inputs  
        bb: BlackBoxIdeal
             an black box ideal
    Outputs
        : HashTable
            with  keys { \tt precision } (the required length of the computed jets)
            and { \tt numTrials } ( the number of required jet computations )
    Description
        Text
            Show currently used parameters for @TO2{isCertainlySingularAt,"point singularity test"}@
            There are two parameters, the required length of the computed jets ( { \tt precision } )
            and  the number of required jet computations ({ \tt numTrials }). If one of the jet computations fails,
            the point is certainly not smooth.

            Let us construct a black box
        Example
            R = QQ[x,y]
            bbI = blackBoxIdeal ideal(x^2-y^3);
        Text
            Now we may look at the default singularity test parameters:
        Example
            bbI.singularityTestOptions()
        Text
            The test parameters may be modified by the user with @TO{setSingularityTestOptions}@
        Example
            jetLength = 3;
            numTrials = 5;
            bbI.setSingularityTestOptions(jetLength, numTrials);
            bbI.singularityTestOptions()
    SeeAlso
        setSingularityTestOptions
        isProbablySmoothAt
        isCertainlySingularAt
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
          trying to find a @TO2{jetAt,"jet"}@ starting at the point.
          
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
          For the tests the length is taken by the precision of 
          @TO{singularityTestOptions}@ and the required number of successful trials  
          is given by 'numTrials'
        Example    
          bbI.singularityTestOptions()
        Text
          The default singularityTestOptions can be changed with
          @TO{setSingularityTestOptions}@
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
          the black box ideal, an exception PointNotOnBlackBox is thrown
        Example
          pointNotOnCurve = matrix{{4,5_QQ}}
          bbI.isZeroAt(pointNotOnCurve)
          catch bbI.isCertainlySingularAt(pointNotOnCurve)
          wantedJetLength = 1
          catch jetAt(bbI,pointNotOnCurve,wantedJetLength)
    SeeAlso
       setSingularityTestOptions
       jetAt
      
///

doc ///
    Key
        jetAtOrException
        (jetAtOrException, BlackBoxParameterSpace, Matrix, ZZ)
        "jetAt"
    Headline
        finds a jet on a variety defined by a black box
    Usage   
        jetAt(bb,point,length)
        jetAtOrException(bb,point,length)
        bb.jetAt(point,length)
        bb.jetAtOrException(point,length)
    Inputs  
        bb: BlackBoxIdeal
             an black box ideal
        point: Matrix 
             the coordinates of a point
        length: ZZ
             the length of the desired jet
    Outputs
        : Jet
    Description
        Text 
          "jetAtOrException" can be abbreviated as "jetAt".
          It tries to find a jet starting at a given point on 
          a variety given by a black box.
          
          If the variety defined by the black box is smooth at
          the given point, arbitray jets can be found by the
          implemented algorithm. If the point is singular, the
          algorithm fails after a finite number of steps. If 
          this happens an exception is raised. This is implemented
          using the throw/catch mechanism.
          
          Consinder for example a cuspidal cubic:
        Example
          Fp = ZZ/101
          R = Fp[x,y]
          I = ideal(x^2-y^3)
          bbI = blackBoxIdeal I;
        Text
          Consider a point on the cuspidal cubic different from the origin:
        Example
          point = matrix{{8,4_Fp}}
        Text
          Check whether the point lies on the cuspidal cubic:
        Example  
          bbI.isZeroAt(point)
        Text
          We now look for a jet:
        Example
          j = jetAt(bbI,point,3)
        Text
          The defining equations of the ideal indeed vanish on the jet:
        Example
          sub(I,j)
        Text
          At the origin the cuspidal cubic is singular. Short jets can be found,
          but not long ones.
        Example
          origin = matrix{{0,0_Fp}}
          catch jetAt(bbI,origin,1)  
          catch jetAt(bbI,origin,2)  
          catch jetAt(bbI,origin,3)  
        Text
          Notice that one has to use the catch/throw mechanism
          to obtain readable error messages. 
    SeeAlso
        jetAtOrNull
        isProbablySmoothAt
        isCertainlySingularAt
///

doc ///
    Key
        jetAtOrNull
        (jetAtOrNull, BlackBoxParameterSpace, Matrix, ZZ)
    Headline
        finds a jet on a variety defined by a black box
    Usage   
        jetAtOrNull(bb,point,length)
        bb.jetAtOrNull(point,length)
    Inputs  
        bb: BlackBoxIdeal
             an black box ideal
        point: Matrix 
             the coordinates of a point
        length: ZZ
             the length of the desired jet
    Outputs
        : Jet
            or null.
    Description
        Text 
          @TO jetAtOrNull @ does the same as @TO jetAtOrException @,
          i.e. it
          tries to find a jet starting at a given point on 
          a variety given by a black box.          
          The only difference is, that if it failes to
          find a jet of the desired length, no exception is raised.
          Instead null is returned.
          
          Consinder for example a cuspidal cubic:
        Example
          Fp = ZZ/101
          R = Fp[x,y]
          I = ideal(x^2-y^3)
          bbI = blackBoxIdeal I;
        Text
          At the origin the cuspidal cubic is singular. Short jets can be found,
          but not long ones.
        Example
          origin = matrix{{0,0_Fp}}
          jetAtOrNull(bbI,origin,1)  
          jetAtOrNull(bbI,origin,2)  
          catch jetAtOrException(bbI,origin,2)  
        Text
          Notice that with @TO jetAtOrNull @ one does
          not need to use the catch/throw mechanism, but
          also one does not get useful error messages.
    SeeAlso
        jetAtOrException
///

doc ///
    Key
        continueJetOrException
    Headline
        increases the length of a given jet on a variety defined by a black box
    Usage   
        continueJetOrException(bb,jet,length)
    Inputs  
        bb: BlackBoxIdeal
             THIS SHOULD NOT BE NECESSARY!
        jet: Jet 
        length: ZZ
             the length of the desired jet
    Outputs
        : Jet
    Description
        Text 
          This takes a jet of a given length and tries
          to continues the algorithm of @TO jetAt@ until
          the jet has the desired length.
          
          If the variety defined by the black box is smooth at
          the given point, arbitray jets can be found by the
          implemented algorithm. If the point is singular, the
          algorithm fails after a finite number of steps. If 
          this happens an exception is raised. This is implemented
          using the throw/catch mechanism.
          
          Consinder for example a cuspidal cubic:
        Example
          Fp = ZZ/101
          R = Fp[x,y]
          I = ideal(x^2-y^3)
          bbI = blackBoxIdeal I;
        Text
          Consider a point on the cuspidal cubic different from the origin:
        Example
          point = matrix{{8,4_Fp}}
        Text
          Check whether the point lies on the cuspidal cubic:
        Example  
          bbI.isZeroAt(point)
        Text
          We now look for a jet:
        Example
          j = jetAt(bbI,point,3)
        Text
          Now we increase the length of the jet
        Example
          continueJetOrException(bbI,j,4)
        Text
          At the origin the cuspidal cubic is singular. Short jets can be found,
          but not long ones.
        Example
          origin = matrix{{0,0_Fp}}
          j = jetAt(bbI,origin,1)  
          catch continueJetOrException(bbI,j,3)
        Text
          Notice that one has to use the catch/throw mechanism
          to obtain readable error messages. 
    SeeAlso
        jetAtOrException
///


doc ///
    Key
        jetStatsAt
    Headline
        counts possible jet lengths at a singular point
    Usage   
        jetStatsAt(bb,point,trials,length)
        bb.jetStatsAt(point,trials,length)
    Inputs  
        bb: BlackBoxIdeal
             a black box ideal
        point: Matrix 
             the coordinates of a point
        trials: ZZ
             the number of times the jet-finding procedure
             is started
        length: ZZ
             the maximum length used in the search
    Outputs
        : Jet
    Description
        Text 
          "jetAtOrException" can be abbreviated as "jetAt".
          It tries to find a jet starting at a given point on 
          a variety given by a black box.
          
          If the variety defined by the black box is smooth at
          the given point, arbitray jets can be found by the
          implemented algorithm. If the point is singular, the
          algorithm fails after a finite number of steps. If 
          this happens an exception is raised. This is implemented
          using the throw/catch mechanism.
          
          Consinder for example a cuspidal cubic:
        Example
          Fp = ZZ/5
          R = Fp[x,y]
          bbCusp = blackBoxIdeal ideal(x^2-y^3);
        Text
          At the origin the cuspidal cubic is singular. Short jets can be found,
          but not long ones.
        Example
          origin = matrix{{0,0_Fp}}
          catch jetAt(bbCusp,origin,3)  
        Text
          Lets check how often this happens:
        Example
          jetStatsAt(bbCusp,origin,10,5^3)
        Text
          Now lets do the same for a node:
        Example
          bbNode = blackBoxIdeal ideal(x*y);
          jetStatsAt(bbNode,origin,10,5^3)
        Text
          Notice that the statistics is significantly different
          form the cusp example. Possibly some help in 
          classifying implicitly given singularities can be
          obtained from this.  
    SeeAlso
        jetAtOrException
        isProbablySmoothAt
        isCertainlySingularAt
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
          that can evaluate the generators of the black box ideal. A black box parameter space
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
           bbC.pointProperties()
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
          bbC.pointProperties()
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
          so that for example pointProperties will return
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
          bbI.pointProperties()
        Text
          A blackBoxParameterSpace usually does not have 
          the property "jacobianAt" since such a parameter space
          does not even have an implicit representation of 
          the equations.
        Example
          bbParam = blackBoxParameterSpace(2,QQ);     
          bbParam.hasPointProperty("jacobianAt") 
          bbParam.pointProperties()
        Text
          To illustrate why this can happen think for example of 
          the space of
          singular cubics. In a BlackBoxParameterSpace one
          would simply test the smoothness of a give cubic
          via Groebner basis calculation. This does not automatically
          give rise to a representation of the corresponding
          Diskriminant. 
    SeeAlso
         pointProperties
         registerPointProperty
///

doc ///
    Key
        "attributes"
        "numVariables"
    Headline
        list the attributes of a black box 
    Usage   
        bbI.attributes()
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
          bbI.attributes()
        Text
          Lets look at these attributes in turn
        Example
          bbI.coefficientRing
          bbI.equations
          bbI.ideal
          bbI.jacobian
          bbI.numVariables
          bbI.ring
          class bbI
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
          bbE.attributes()
        Text
          Notice that "equations" and "jacobian" are missing, since
          no explicit equations of the blackBoxIdeal are provided.
        Example
          bbP = blackBoxParameterSpace(2,QQ);
          bbP.attributes()
        Text
          For a blackPointParameterSpace "ring" and "unknowns" are 
          missing since there are no equations (not even implicit ones)  
    SeeAlso
         pointProperties
         memberMethods
///

doc ///
    Key
        "memberMethods"
        (memberMethods, BlackBoxParameterSpace)
    Headline
        list the methods of a black box 
    Usage   
        bbI.memberMethods()
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
          bbI.memberMethods()
        Text
          Each of these methods has their own documentation node.
    SeeAlso
         hasPointProperty
         attributes
         pointProperties         
         registerPointProperty
         rpp
         setSingularityTestOptions
         updatePointProperty
         upp
///

doc ///
    Key
        "pointProperties"   
        (pointProperties, BlackBoxParameterSpace)
    Headline
        list the pointProperties of a black box 
    Usage   
        bbI.pointProperties()
        pointProperties bbI
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
          bbI.pointProperties()
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

doc ///
    Key
        "setSingularityTestOptions"
    Headline
        change how singularities are detected
    Usage   
        bbI.setSingularityTestOptions(precision,numTrials)
    Inputs  
        bbI: BlackBoxIdeal
        precision: ZZ
             length of jets used
        numTrials: ZZ
             number of such jets required     
    Outputs
        : Boolean
    Description
        Text
          To test wether a given point P on a variety given by
          a black box ideal is smooth, the package tries to construct
          several jets of certain length. The number of required trials and jetlengths
          is given by @TO{singularityTestOptions}@ and the default is to compute
          2 jets of length 10 starting at P. If the
          variety is smooth at P such jets always exists. If the variety is
          singular at P a generic jet can not be extended to arbitrary length.
          If the required number of jets can not be found the point property 
          isCertainlySingular has the value true. If the required number 
          of jets can be found the point property @TO{isProbablySmoothAt}@ has the
          value true.
        Example
          K = QQ
          R = QQ[x,y,z]      
          I = ideal (x*z,y*z)
          bbI = blackBoxIdeal I;
          pointOnLine = matrix{{0,1,0_K}}
          bbI.isProbablySmoothAt(pointOnLine)
          origin = matrix{{0,0, 0_K}}
          bbI.isProbablySmoothAt(origin)
          bbI.isCertainlySingularAt(origin)
        Text
          Since the singularity of the above curve
          at the origin is of small degree
          it is detected by looking at jets of length 10.
          If we change the test options to look only at jets of length 4,
          the singularity can not be detected.
        Example
          bbI.setSingularityTestOptions(4,1)
          bbI.singularityTestOptions()
          bbI.isCertainlySingularAt(origin)
        Text
          The construction of jets is time intensive. For many applications
          precision=2 and numTrials=1 is sufficient, even if this only detects
          the most simple singularities.  
    Caveat
          This works only with black box ideals, since they contain an algorithm
          that can calculate the derivative of the equations at a given point
          (jacobianAt). This is needed to construct jets iteratively.
    SeeAlso
          singularityTestOptions
          isCertainlySingularAt
          isProbablySmoothAt
          jetAt
///

doc ///
    Key
        "Jet"
    Headline
        a type for handling jets
    Description
       Text
          Algebraically a jet is a map 
          
          R \to K[e]/e^{d+1} 
          
          where R is a polynomial ring in n variables.
          Geometrically it is a truncated curve germ.
          
          Here a jet is modeled as a row matrix with n
          entries in an eps-Ring of precision d (see @TO getEpsRing @).
          
          As an example consider the cuspidal cubic curve
          in the plane:                    
       Example
          K = ZZ/7
          R = K[x,y]      
          I = ideal (x^2-y^3)
          bbI = blackBoxIdeal I;
       Text
          the following point lies on the cuspidal cubic
       Example
          smoothPoint = matrix{{1,1_K}}
          bbI.isZeroAt(smoothPoint)
       Text
          Lets now make a jet of lenght 2 lying on cuspidal
          cubic and starting in the above point.
       Example   
          j = bbI.jetAt(smoothPoint,2)
       Text
          The jet does indeed lie on the cuspidal cubic
       Example
          sub(I,j)
       Text
          A jet remembers from which point it started
          from:
       Example
          j#"point"    
       Text
          Also it knows its length:
       Example
          j#"jetLength"
       Text
          Finally it remembers from which BlackBoxIdeal it
          was created
       Example
          class j#"parent"
       Text
          The length of a jet can be increased:
       Example
          j3 = continueJet(bbI,j,3)    
       Text
          Jets can only be reliably created in smooth
          points of a variety (by a variant of Newtons method).
          If the starting point of the jet is not smooth
          on the variety defined by the BlackBoxIdeal, then
          after a finite number of steps the algorithm can 
          not continue. If this happens
          an exception is raised. 
          
          For readable error messages
          and more precise information one has to use 
          the throw/catch mechanism:
       Example
          singularPoint = matrix{{0,0_K}};
          catch bbI.jetAt(singularPoint,1)          
          catch bbI.jetAt(singularPoint,2)
       Text                 
          Notice that the search for fails at length 2 most of the time,
          since the singularity has multiplicity 2. If one tries
          long enough, a longer jet can be found (lying on one
          of the branches at the origin):
       Example        
          jetStatsAt(bbI,singularPoint,3,200) 
    Caveat
    SeeAlso
///

doc ///
    Key
        "JetSet"
    Headline
        a type for handling a set of jets starting at the same point
    Description
       Text
          A JetSet is a set of jets that start at the same point.
          This is important since jets starting at the
          same (smooth) point of a variety X must 
          lie on the same component of X.
       
          For example consider the cuspidal cubic:
       Example
          K = ZZ/101
          R = K[x,y]
          bbCusp = blackBoxIdeal ideal(x^2-y^3);
       Text
          Lets make a jet at a smooth point
       Example
          smoothPoint = matrix{{1,1_K}}
          j = jetAt(bbCusp,smoothPoint,3)
       Text
          We now make a JetSet containing one Element
       Example
          js = new JetSet from j
       Text
          We add another jet to this collection
       Example
          addElement(js,jetAt(bbCusp,smoothPoint,3)) 
       Text
          Lets now consider a differnent point on 
          the cuspidal cubic
       Example
          otherSmoothPoint = matrix{{8,4_K}}
          otherJet = jetAt(bbCusp,otherSmoothPoint,3)
       Text
          this can not be added to our JetSet, because
          the jet starts at a different point.

          -- addElement(js,otherJet)
    Caveat
    SeeAlso
///




TEST ///
    --test for issue #117
    coeffRing := ZZ/3;
    numVariables := 2;
    R = coeffRing[x,y];

    bb = blackBoxParameterSpace( numVariables , coeffRing );

    P = matrix{{0_coeffRing, 0_coeffRing}}
    valuesAt := (blackBox, point)-> matrix {{0_coeffRing }};
    bb = bb.rpp("valuesAt",valuesAt);

    bb.valuesAt(P)
    bb.isZeroAt(P)
    jetAt(bb,P,1)
///

TEST ///
    -- test for issue #111
    --
    coeffRing := ZZ/3;
    numVariables := 2;
    bb = blackBoxParameterSpace( numVariables , coeffRing );

    try (bb.setSingularityTestOptions(5,5);) then (error "should fail";) else();

    try (bb.singularityTestOptions();) then (error "should fail";) else();

    valuesAt := (blackBox, point)-> matrix {{5}};
    bb = bb.rpp("valuesAt",valuesAt);
    bb.singularityTestOptions()
    bb.setSingularityTestOptions(5,5);
    point = matrix{{0_coeffRing, 1_coeffRing}}
    bb.valuesAt (point)
    catchOrResult = catch bb.isCertainlySingularAt(point)
    assert (class catchOrResult === PointNotOnBlackBox)
///

TEST ///
    --test for issue #86
    kk = ZZ/2
    R = kk[x,y]
    I = ideal(x^5-y^7);
    bbI = blackBoxIdeal I;
    smoothPoint = matrix{{8,4_kk}}
    bbI.valuesAt(smoothPoint)
        
    bbI.isProbablySmoothAt(smoothPoint)
    origin = matrix{{0,0_kk}}
    bbI.isProbablySmoothAt(origin)
    bbI.isCertainlySingularAt(origin)
    bbI.setSingularityTestOptions(4,1)
    bbI.isCertainlySingularAt(origin)
    bbI.isProbablySmoothAt(origin)
///


TEST ///
    -- test issue #21
    
    coeffRing := ZZ/3;
    numVariables := 2;
    bb = blackBoxParameterSpace( numVariables , coeffRing );

    R = coeffRing[x,y];
    valuesAt := (bb, point)-> matrix {{sub(  x^2, point )}};
    bb = bb.rpp("valuesAt",valuesAt);
    P = matrix{{0_coeffRing, 0_coeffRing}}
    bb.jacobianAt(P)
    P = matrix{{1_coeffRing, 1_coeffRing}}
    valuesAt(bb,P);
    
    --bb.jacobianAt(point) should throw an  PointNotOnBlackBox exception
    assert (PointNotOnBlackBox === class catch  bb.jacobianAt(P) )     

    valuesAt = ( point)-> matrix {{sub(  x^2, point )}};
    valuesAt(P);
    bbE = blackBoxIdealFromEvaluation (R , valuesAt );

    P = matrix{{0_coeffRing, 0_coeffRing}}
    bbE.jacobianAt(P)
    P = matrix{{1_coeffRing, 1_coeffRing}}
    --bbE.jacobianAt(point) should throw an  PointNotOnBlackBox exception
    assert (PointNotOnBlackBox === class catch  bbE.jacobianAt(P) ) 
    
    
    C = QQ;
    R = QQ[x];
    bbI = new BlackBoxIdeal from ideal x;
    point = matrix{{1_QQ}};
    
    -- bbI.jacobianAt(point) should gove an error , becaues the point is not on the variety
    assert(PointNotOnBlackBox === class catch  bbI.jacobianAt(point) )      
    

///


TEST ///
    -- test issue #35
    K = QQ;
    R = K[x,y];
    bbI = new BlackBoxIdeal from ideal x*y;
    point = matrix{{1_QQ, 1_QQ}};

    -- not on black box throws error
    jetOrError = catch jetAt(bbI,point,1)
    assert (class jetOrError === PointNotOnBlackBox)
    
    point = matrix{{0_QQ, 0_QQ}};
    jet =  jetAt(bbI,point,1)    
    assert(jet =!= null)
    assert(length jet == 1)
    jet =  jetAt(bbI,point,0)
    assert(length jet == 0)
///



-- JK we have to put the undocumented statement at the end, because undocumented checks, if the 
-- functions or symbols are already defined!
--
-- the following symbols which are marked as undocumented are in fact documented 
-- inside the BlackBoxParameterSpace and BlackBoxIdeal
-- please do only mark documented symbols as undocumented, 
-- at least there should be a comment note inside this package.
-- 
undocumented { 
    (isCertainlySingularAt, BlackBoxParameterSpace,Matrix,HashTable),
    (isCertainlySingularAt, BlackBoxParameterSpace, Matrix, ZZ, ZZ),
    setJetLengthHeuristic,
    setMonomialDegreeHeuristic,
    setName,   
    enableChecks,
    disableChecks,
    withChecks,
    (net,Jet),
    (net,JetSet),
    (net,InterpolatedIdeal),
    (net, BlackBoxParameterSpace),  
    setIsZeroAt,      -- internal   
    setJacobianAt,    -- internal   
    setValuesAt,      -- internal   
    setPointProperty, -- internal   
    --(deduceNumGenerators), -- internal   
    dropDegreeInfo,      -- internal   
    updateBlackBox,  --internal
    equations,
    keysWithoutSymbols,
    checkInputPoint,
    pointPropertiesAsSymbols,
    unknownIsValid,
    unknowns,
    numTrials,
    pointProperty, -- internal function
    guessAcceptedParameterNumber, -- internal function
    updateSingularityTest, --internal function
    createMapHelper,
    JetLengthHeuristic,
    ConstantJetLengthHeuristic,
    constantJetLengthHeuristic,
    BasicJetLengthHeuristic,
    basicJetLengthHeuristic,
    InterpolationMonomialDegreeHeuristic,
    BasicInterpolationMonomialDegreeHeuristic,
    ConstantInterpolationMonomialDegreeHeuristic,
    setMonomialDegreeHeristic,
    interpolationTargetJetLength,
    setJetLength,
    setAdditionalJetLength,   
    targetMonomialDegree,
    clearCache,
    minComponentDegree,
    maxInterpolationDegree,
    maxComponentDegree,
    blackBox,
    increaseInterpolationJetLength,
    decreaseInterpolationJetLength,
    setSameComponentPrecision,
    setComponentNamePrefix,
    sameComponentTargetJetLength,
    additionalJetLength,
    interpolator,
    componentNameInUse,
    componentNamesInUse,
    jetSet,
    
    maximalConditions,    
    "InterpolatedIdeal ? InterpolatedIdeal", -- implemented comparison by name for sorting purposes. probably not a good idea, since the order holds for the names but not for the ideals !!
    -- purpose of overriding new from Thing is to disallow arbitrary HashTables as objects
    (NewFromMethod, JetSet, Thing),
    (NewFromMethod, BlackBoxIdeal, Thing),
    (NewFromMethod, BlackBoxParameterSpace, Thing),
    (NewFromMethod, JetLengthHeuristic, Thing),
    (NewFromMethod, BasicInterpolationMonomialDegreeHeuristic, Thing),
    (NewFromMethod, BasicJetLengthHeuristic, Thing),
    (NewFromMethod, ConstantInterpolationMonomialDegreeHeuristic, Thing),
    (NewFromMethod, ConstantJetLengthHeuristic, Thing),
    (NewFromMethod, InterpolatedIdeal, Thing),
    (NewFromMethod, InterpolationMonomialDegreeHeuristic, Thing),
    (NewFromMethod, MapHelper, Matrix),
    transformedAnswer
} 



end

-- undocumented: "compatible"


uninstallPackage"BlackBoxIdeals"
loadPackage "BlackBoxIdeals"
installPackage"BlackBoxIdeals"
check BlackBoxIdeals




-- Todo: introduce force mode option for registerPointProperty? (overwriting registered property?)
-- todo: introduce for all known properties a precomposition with checkInputPoint
    


-- Estring = new Type of String ; the goal was to distinguies string keys 
-- from symbol keys without changing string visuaulization (net())
-- but that turned out to be impossible(or too hard), because e.g. the operator '#' cannot be overloaded.


--net (String) := Net =>(str)->(   strl := "\""| str| "\"" ;   return stack separate strl; );


--(JK) why did I need comparison between a string and a symbol ? (jk) - if I recall correctly, it was related for intervals...

String ? Symbol := (str,symb)->
(
    1 ? 2
);

Symbol ? String := (str,symb)->
(
    2 ? 1
);
-- Q: how to get the parameter list of a function or method?
--
-- answer:  not possible yet, but it is possible to get the number of parameters 
-- using 'disassemble' for functions (where the number of variables is fixed and not variable)
-- and for methods it forks as follows: consider a method foo with several functions installed.
-- then   apply(methods foo, m-> (m, disassemble lookup m) )
   
-- just after exporting:
-- why is/was setting 'jacobianAt' necessary??? (jk)

-- jacobianAt := global jacobianAt;
-- jacobianAt = global jacobianAt;
--
--

-- jetAtWrapper currently deprecated !
-- options makes more sense for trials (=1)
jetAtWrapper = method( Options=>{"jetLength" =>4, "trials"=>5} );

jetAtWrapper( BlackBoxParameterSpace, Matrix )  := MutableHashTable => opts -> ( blackBox,  point )  ->
(
    return jetAt( blackBox,  point, opts#"jetLenth", opts#"trials");
);



-- need test for randomIterator (for fixed error : to less trials  )
--@TO2{jetAt,"jet"}@ computations at a point independently of the ideal representation. \break \break 


--        (clearCoeffDenominators, Ideal )    


-- dient dazu, an einen Punkt die BlackBox mit anzuheften.
Point = new Type of HashTable;

pointObject = method();
pointObject (Thing, Matrix) := Point =>(parent, point)->
(
    resultPoint := new HashTable from {("parent", parent),
                                       ("point", point)};
    resultPoint 
)


-- parent is needed for coercion

Point := new Type of HashTable;


pointObject (Thing, Matrix) := (thing, point)->
(
    if not (thing.hasElement(point)) then 
        error(" point is not on parent ");
    pointObject = new HashTable from {("parent",thing),
                                      ("point",point)
                                     }
    return pointObject;
);


BlackBoxPoint := new Type of Point
{
    parent is of type BlackBoxParameterSpace or derived
    point is a Matrix?
}

blackBoxPointObject (BlackBoxParameterSpace, Matrix) := (blackBox, point)->
(
    if not (blackBox.isZeroAt(point)) then error(" point is not on BlackBox ");
    pointObject = new HashTable from {("parent",blackBox),
                                      ("point",point)
                                     }
    return pointObject;
);

-- single trial returns 

(jet, bestJet, failedJetLength)





JetInfoAt 
(
    bestJet
    worstJet
    Tally failedJetLength=> count
    targetLength
)

-- or even better, JetInfoAt contains pairs length -> list of jets with that length.

-- 
 
-- and now we can pass the black box and the point ! (we have free another two parameters)

-- we could do precision Jet and then check  first "jetLength", then 'precision ring jet'.


InterpolatedIdeal: parent is JetSet

Jet: parent is Point ? (or black box and point)




worstJetAt = method();

worstJetAt( BlackBoxParameterSpace, Matrix, ZZ, ZZ) := Jet => ( blackBox,  point, jetLength, numTrials )  ->
(
    if ( numTrials<1 ) then error "jetAt: expected numTrials >=1 ";
    
    jetResult  := jetAtSingleTrial ( blackBox,  point, jetLength);      
    
    worstJet := jetResult#"bestJet";
    
    for i in 2..numTrials do
    (
        newJetResult := jetAtSingleTrial ( blackBox,  point, jetLength);

        if ( length worstJet > length newJetResult#"bestJet" ) then
        ( 
            worstJet = newJetResult#"bestJet";
        );       
    );
    return worstJet;
);


-- returns the longest jet with length <=jetLength obtained after maximal numTrial trials
bestJetAt= method();

bestJetAt( BlackBoxParameterSpace, Matrix, ZZ, ZZ) := Jet  => ( blackBox,  point, jetLength, numTrials )  ->
(
    if ( numTrials<1 ) then error "jetAt: expected numTrials >=1 ";
    
    jetResult := jetAtSingleTrial ( blackBox,  point, jetLength);    

    if  jetResult#"succeeded"  then  return jetResult#"jet"; 
    
    bestJet := jetResult#"bestJet";

    for i in 2..numTrials do
    (
        newJetResult  := jetAtSingleTrial ( blackBox,  point, jetLength);

        if  newJetResult#"succeeded"  then return  newJetResult#"jet"; 

        if ( length bestJet  < length  newJetResult#"bestJet"  ) then
        ( 
            bestJet =  newJetResult#"bestJet";
        );       
    );
    return bestJet;
);
--locus  (JetSet) := Thing => (jetSet )->
--(
--    return locus firstElem jetSet;
--);
