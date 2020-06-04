

bblog := BlackBoxLogger;
-- InterpolatedComponent = new Type of MutableHashTable;

BlackBoxInterpolator = new Type of HashTable;


OnComponentAnswerStrategy = new Type of HashTable;




NullIfNotSmooth = new Type of OnComponentAnswerStrategy;
ExceptionIfNotSmooth = new Type of OnComponentAnswerStrategy;
SmoothnessInfoWithAnswerPair = new Type of OnComponentAnswerStrategy;
PlainTextSmoothnessInfoWithAnswerPair = new Type of OnComponentAnswerStrategy;

nullIfNotSmoothStrategy = ()->
(
    answerTransformer := new MutableHashTable;
    answerTransformer.transformedAnswer = (isSmooth, answer) ->
    (
        if (not isSmooth) then return null;
        return answer;
    );
    answerTransformer = newClass(NullIfNotSmooth, answerTransformer);
    return answerTransformer;
);

exceptionIfNotSmooth = ()->
(
    answerTransformer := new MutableHashTable;
    answerTransformer.transformedAnswer = (isSmooth, answer) ->
    (
        if (not isSmooth) then throw new SingularPointException;
        return answer;
    );
    answerTransformer = newClass(ExceptionIfNotSmooth, answerTransformer);
    return answerTransformer;
);

smoothnessInfoWithAnswerPair = ()->
(
    answerTransformer := new MutableHashTable;
    answerTransformer.transformedAnswer = (isSmooth, answer) ->
    (        
        return (isSmooth, answer);
    );
    answerTransformer = newClass(SmoothnessInfoWithAnswerPair, answerTransformer);
    return answerTransformer;
);

smoothPointText := "probablySmooth";
singularPointText := "certainlySingular";

plainTextSmoothnessInfoWithAnswerPair = ()->
(
    answerTransformer := new MutableHashTable;
    answerTransformer.transformedAnswer = (isSmooth, answer) ->
    (       
        if isSmooth then 
            return (smoothPointText, answer);
        return (singularPointText, answer);
    );
    answerTransformer = newClass(PlainTextSmoothnessInfoWithAnswerPair, answerTransformer);
    return answerTransformer;
);






eassert = method();

eassert(Boolean,String) := Nothing =>(statement, errorMessage)->
(
    if (not statement) then 
        error errorMessage;
);

monomialBasisSize = method();

monomialBasisSize (ZZ, ZZ , Ring) := ZZ => (minDegree, maxDegree, ring)->
(
    eassert(maxDegree>=0, "maxDegree was smaller 0");
    eassert(minDegree>=0, "minDegree was smaller 0");
    eassert(minDegree<=maxDegree, "maxDegree was smaller minDegree");
    degrees := apply(maxDegree-minDegree+1, i->minDegree+i);
    basisSize := sum apply(degrees, deg->rank source super basis(deg,ring) );
    return basisSize;
);

monomialBasisSize (ZZ , Ring) := ZZ => ( maxDegree,ring)->
(
     eassert(maxDegree>=0, "maxDegree was smaller 0");
    minDegree := 0;
    return monomialBasisSize(minDegree, maxDegree,ring);
);



TEST ///
    R = ZZ[x,y]
    monomialMaxDegree = 4;
    mons = matrix {flatten apply(monomialMaxDegree+1, currentDegree->flatten entries basis(currentDegree,R))};
    assert (rank source mons == monomialBasisSize(4,R));
   
///

InterpolatedComponent = new Type of HashTable;

new InterpolatedComponent from Thing :=  (InterpolatedIdealAncestor,l)->
(  
    ---type check?
    error("InterpolatedComponent from Thing not supported or not implemented");
);

InterpolatedComponent ? InterpolatedComponent := (i1,i2)->
(
    return (i2#"name"() ? i2#"name"());
)
--new InterpolatedComponent from MutableHashTable :=  (InterpolatedIdealAncestor,l)->
--(  
    ----type check?
    --return l;
--);

--new InterpolatedComponent from List :=  (InterpolatedIdealAncestor,l)->
--(  
--    return new InterpolatedComponent from new MutableHashTable from l;
--);


netInterpolatedComponent = (II)->
(
    col0  := {"\"" | II#"name"() | "\","};
     
    col1  := {"maxDegree"
              ,"point"
              ,"jetSet" 
              --,"blackBox"
             };
  
    col2 :=  {"=>"
              ,"=>"
              ,"=>"
              --,":"
              };
    
    col3  :=  {
             net II#"maxDegree"
             ,net II#"point"
             ,net class II#"jetSet" | "{.."| toString size II#"jetSet" |"..}"
             --,net class II#"blackBox"
             };

    stackCol0 := stack col0;
    stackCol1 := stack col1;
    stackCol2 := stack col2;
    stackCol3 := stack col3;
    
    stackBottomRight :=  stackCol1 | " " |stackCol2 | " " |stackCol3;
    
    stackRight := stack { net II#"ideal" , stackBottomRight };
    
    result := stackCol0 | " " | stackRight;
    bracedResult:= net class II | "{" | result | "}";
    return bracedResult;
);

net (InterpolatedComponent) := Net =>(ii)->
(
    return netInterpolatedComponent(ii);
)


-- find polynomials containing a component
-- via interpolation

-- interpolate()
--
-- find linear combinations of monomials containing a given jet
-- it seems that currently jetList have only one entry.
--

interpolate = method();

interpolate (Matrix, List) := Ideal => (mons, jetList) -> 
(
    bblog.debug("--interpolate call;");
    bblog.debug("--interpolate call; \n --mons : " | toString mons );
    bblog.debug("--interpolate call;--jetList :" | toString jetList);
    bblog.debug("--interpolate call ok");    
    R := ring mons;
    K := coefficientRing R;
    failedJets := select(jetList, jetP->jetP===null);    
    if #failedJets >0 then 
    ( 
        error throw new SingularPointException from {"errorMessage"=>"jets not succeed. Point is not smooth?";}
    );
    --
    -- substitute jets into the monomials and take coefficients
    -- die jets werden in jedes monom eingesetzt. (matrix : #monome x 1 für einen jet.)
    -- substitutedList := apply(jetList, jetP ->   sub(mons, jetP#"value"));
    -- nun haben wir für je einen jet J anstatt n monome n polynome in eps. 
    -- 'last coefficients' gruppiert die Koeffizienten nach Grad von 'eps', d.h. für jeden Grad von eps 
    -- haben wir eine Matrix-Zeile mit #monome Einträgen.
    -- BlackBoxLogger.debug("interpolate, substituded:" | toString substitutedList);
    -- coeffs := apply(substitutedList, sp-> coefficients sp);    
    -- BlackBoxLogger.debug("interpolate, coeffs:" | toString coeffs);
    --
    coeffList := apply(jetList,jetP ->  sub(last coefficients sub(mons, jetP#"value"),K));
    BlackBoxLogger.debug("interpolate, coeffList:" | toString coeffList);    
    --
    -- haben wir einen weiteren Jet, so muss die gleiche Linearkombination der Spalten für den ersten
    -- auch für den zweiten Jet gelten
    --  => wir hängen die Koeffizientenmatrix für den zweiten eingesetzten jet einfach an die erste Matrix unten dran
    --    (weitere Zeilen kommen hinzu, Kommando 'fold'; mit || als Konkatenation)
    --
    folded := fold((a,b)->(a||b),coeffList);
    BlackBoxLogger.debug("interpolate, folded coefflist" | toString folded);
    -- nun suchen wir für die Spalten eine Linearkombination, so dass jede Zeile zu 0 wird. (jet auf der Komponente).
    -- find interpolation solution    
    s := syz folded; -- todo question: is the syz command necessary here? 
    BlackBoxLogger.debug("interpolate, syzygies" | toString s);
    -- make polynomials from the solution
    -- (jk) todo: is mingens fast enough?
    I := ideal mingens ideal(mons*s);
    return I;
)

 

createInterpolatedComponent = method();
createInterpolatedComponent( Ideal, ZZ, Thing, BlackBoxParameterSpace ) := InterpolatedComponent => 
                       (I,maxDegree, jetset, blackBox)->
(

    point := jetset#"point";
    
    name := "undefined";
    
    getName := ()->
    (      
        return name;
        --return blackBox.componentNameByBasePoint(point);
    );
    
    -- chec for component name conflicts:
    if (  blackBox.componentNameInUse(name)) then
    (
         error ("name '" | name | "' is already used by a differen component");
    ); 
    
    setName := (newName)->
    (
       --get components, check if there is no conflict, then rename.     
        if (newName == name ) then return ; --nothing changed
      
        -- check for component name conflicts:
        if (  blackBox.componentNameInUse(newName)) then
        (
            error ("name '" | newName | "' is already used by a differen component");
        )       
        else
        (
            name = newName;
            -- when we here notify an observer 
            -- the  observer could  maintain a cache table for (component,name) pairs
        );        
    );   
       
    
     result := new HashTable from {    
        "blackBox" => blackBox,        
        "ideal" => I,        
        "maxDegree" => maxDegree,
        "jetSet" => jetset,
        "point" => point,
        "name" => getName,
        "setName" => setName,
        --(global blackBox) => blackBox,
        --(global ideal) => I,
        --(global maxDegree) => maxDegree,
        --(global jetSet) => jetSet,
        --(global point) => point,
        (global name) => getName,
        (global setName) => setName,
        
    };
    result = newClass(InterpolatedComponent, result);
    return result;
);

createInterpolatedComponent( Ideal, ZZ, String ) := InterpolatedComponent => 
                       (I,maxDegree,  description)->
(
    return createInterpolatedComponent (I,maxDegree, null, description);
);

ideal ( InterpolatedComponent) := Ideal => (ii)->
(
    return ii#"ideal";
)

 
TEST ///
    -- loadPackage "BlackBoxIdeals"
    R = ZZ[x,y]
    I = ideal (x*y)
    bb = new BlackBoxIdeal from I
    p = matrix {{1,0_QQ}}
    jetset = new JetSet from bb.jetAt(p,1)
    II = createInterpolatedComponent(I, 1, jetset, bb)
    ideal II 
///

 



-- createSimpleInterpolator should be private, since each blackbox needs an own component calculator.

-- or, the jets should be stored in the blackBox


JetLengthHeuristic = new Type of HashTable;

new JetLengthHeuristic from Thing := (E,thing)->
(
    error "not supported";
);


ConstantJetLengthHeuristic = new Type of JetLengthHeuristic;

new ConstantJetLengthHeuristic from Thing := (E,thing)->
(
    error "not supported";
);
constantJetLengthHeuristic = method();
constantJetLengthHeuristic (BlackBoxInterpolator, ZZ):= ConstantJetLengthHeuristic => (interpolatorP, targetJetLength)->
(
    interpolator := interpolatorP;
    
    jetLengthEstimator:= new MutableHashTable;
    
    jetLengthEstimator.setJetLength = (newTargetJetLength)->
    (
        targetJetLength = newTargetJetLength;
    );
        
        
    jetLengthEstimator.interpolationTargetJetLength = method();
     
    -- for interpolating at a point
    --
    jetLengthEstimator.interpolationTargetJetLength ( Matrix, ZZ ) :=  ZZ => 
        (   point, monomialDegree )->
    (
        return targetJetLength;
    );
    
    -- probably we do not need this one
    jetLengthEstimator.interpolationTargetJetLength (InterpolatedComponent, Matrix ) :=  ZZ => 
        ( interpolatedIdeal, point )->
    (
        return targetJetLength;
    );
    
     jetLengthEstimator.interpolationTargetJetLength ( InterpolatedComponent,  ZZ  ) :=  ZZ =>
        ( interpolatedIdeal, monomialDegree )->
    (
        return targetJetLength;
    );
    
    jetLengthEstimator.interpolationTargetJetLength (  ZZ  ) :=  ZZ =>( monomialDegree )->
    (
        return targetJetLength;
    );
    
    
    jetLengthEstimator.sameComponentTargetJetLength = method();
    -- component, point (check if point on component)
    --
    jetLengthEstimator.sameComponentTargetJetLength (InterpolatedComponent, Matrix)  :=  ZZ => 
        ( component, point)->
    (
        return targetJetLength;
    );
    
     jetLengthEstimator.sameComponentTargetJetLength (InterpolatedComponent, Matrix, ZZ)  :=  ZZ => 
        ( component, point, maxMonomialDegreeP)->
    (
        return max(targetJetLength,maxMonomialDegreeP);
    );
    
    result := newClass(ConstantJetLengthHeuristic, jetLengthEstimator);
    return result;
);



BasicJetLengthHeuristic = new Type of JetLengthHeuristic;

new BasicJetLengthHeuristic from Thing := (E,thing)->
(
    error "not supported";
);

basicJetLengthHeuristic = method();
basicJetLengthHeuristic (BlackBoxInterpolator) := BasicJetLengthHeuristic => (interpolatorP)->
(
    interpolator := interpolatorP;
    blackBox := interpolator.blackBox;
    
    jetLengthEstimator := new MutableHashTable;
    
    localAdditionalJetLength := 10;
    
    
    
    jetLengthEstimator.setAdditionalJetLength = (additionalJetLengthParam)->
    (
        localAdditionalJetLength = additionalJetLengthParam;
    );
    
    jetLengthEstimator.additionalJetLength = ()->
    (
        return localAdditionalJetLength;
    );
        
          
    jetLengthEstimator.interpolationTargetJetLength = method();
    
     jetLengthEstimator.interpolationTargetJetLength (ZZ) :=  ZZ =>( maxMonomialDegree)->
    (
        basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
        return  basisSize + localAdditionalJetLength;
    );
    
     
    jetLengthEstimator.interpolationTargetJetLength (Matrix,ZZ) :=  ZZ =>( point, maxMonomialDegree)->
    ( 
          basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
          return basisSize + localAdditionalJetLength;
    );
    
    jetLengthEstimator.interpolationTargetJetLength (InterpolatedComponent, Matrix) :=  ZZ =>
         ( interpolatedIdeal, point)->
    (
         maxMonomialDegree := interpolatedIdeal#"maxDegree";        
         basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
         return basisSize + localAdditionalJetLength;
    );
    
    jetLengthEstimator.interpolationTargetJetLength ( InterpolatedComponent, ZZ) :=  ZZ =>
         (interpolatedIdeal, maxMonomialDegree)->
    (
        basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
        return basisSize + localAdditionalJetLength;
    );
    
    
    jetLengthEstimator.sameComponentTargetJetLength = method();
    
     -- component, point (check if point on component)
    jetLengthEstimator.sameComponentTargetJetLength ( InterpolatedComponent, Matrix) :=  ZZ => 
        ( component, point)->
    (
        --maxMonomialDegree := component#"maxDegree";        
        --basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
        --return basisSize + localAdditionalJetLength;
        return interpolator.onComponentPrecision();
    );
    
     jetLengthEstimator.sameComponentTargetJetLength ( InterpolatedComponent, Matrix,ZZ) :=  ZZ => 
        ( component, point, maxMonomialDegreeP)->
    (
        --maxMonomialDegree := max (maxMonomialDegreeP, component#"maxDegree");        
        --basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
        --return basisSize + localAdditionalJetLength;
        return interpolator.onComponentPrecision();
    );
    
    result := newClass(BasicJetLengthHeuristic, jetLengthEstimator);
    return result;
);


jetLengthHeuristicJK = method();
jetLengthHeuristicJK (BlackBoxInterpolator) := BasicJetLengthHeuristic => (interpolatorP)->
(
    interpolator := interpolatorP;
    blackBox := interpolator.blackBox;
    
    jetLengthEstimator := new MutableHashTable;
    
    localAdditionalJetLength := 10;
    
    
    
    jetLengthEstimator.setAdditionalJetLength = (additionalJetLengthParam)->
    (
        localAdditionalJetLength = additionalJetLengthParam;
    );
    
    jetLengthEstimator.additionalJetLength = ()->
    (
        return localAdditionalJetLength;
    );
        
          
    jetLengthEstimator.interpolationTargetJetLength = method();
    
     jetLengthEstimator.interpolationTargetJetLength (ZZ) :=  ZZ =>( maxMonomialDegree)->
    (
        basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
        return  basisSize + localAdditionalJetLength;
    );
    
     
    jetLengthEstimator.interpolationTargetJetLength (Matrix,ZZ) :=  ZZ =>( point, maxMonomialDegree)->
    ( 
          basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
          return basisSize + localAdditionalJetLength;
    );
    
    jetLengthEstimator.interpolationTargetJetLength (InterpolatedComponent, Matrix) :=  ZZ =>
         ( interpolatedIdeal, point)->
    (
         maxMonomialDegree := interpolatedIdeal#"maxDegree";        
         basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
         return basisSize + localAdditionalJetLength;
    );
    
    jetLengthEstimator.interpolationTargetJetLength ( InterpolatedComponent, ZZ) :=  ZZ =>
         (interpolatedIdeal, maxMonomialDegree)->
    (
        basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
        return basisSize + localAdditionalJetLength;
    );
    
    
    jetLengthEstimator.sameComponentTargetJetLength = method();
    
     -- component, point (check if point on component)
    jetLengthEstimator.sameComponentTargetJetLength ( InterpolatedComponent, Matrix) :=  ZZ => 
        ( component, point)->
    (
        maxMonomialDegree := component#"maxDegree";        
        basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
        return basisSize + localAdditionalJetLength;       
    );
    
     jetLengthEstimator.sameComponentTargetJetLength ( InterpolatedComponent, Matrix,ZZ) :=  ZZ => 
        ( component, point, maxMonomialDegreeP)->
    (
        maxMonomialDegree := max (maxMonomialDegreeP, component#"maxDegree");        
        basisSize := monomialBasisSize(maxMonomialDegree,blackBox.ring);
        return basisSize + localAdditionalJetLength;

    );
    
    result := newClass(BasicJetLengthHeuristic, jetLengthEstimator);
    return result;
);


InterpolationMonomialDegreeHeuristic = new Type of HashTable;

new InterpolationMonomialDegreeHeuristic from Thing := (E,thing)->
(
    error "not supported";
);

BasicInterpolationMonomialDegreeHeuristic = new Type of InterpolationMonomialDegreeHeuristic;

new BasicInterpolationMonomialDegreeHeuristic from Thing := (E,thing)->
(
    error "not supported";
);

basicInterpolationMonomialDegreeHeuristic = method();

basicInterpolationMonomialDegreeHeuristic (BlackBoxInterpolator) := 
    BasicInterpolationMonomialDegreeHeuristic => (interpolatorP)->
(
    interpolator := interpolatorP;
    monomialDegreeHeristic := new MutableHashTable;
    
    monomialDegreeHeristic.targetMonomialDegree = method();
    
    
    monomialDegreeHeristic.targetMonomialDegree (Matrix ) := ZZ => (point )->
    (        
        currentMaxDegree := 1;
        
        components := interpolator.componentsAt(point);
        if (#components>0) then 
        (
               currentMaxDegree = max (0, apply (components, c-> max degrees flatten ideal c)) ;  
        );               
        return currentMaxDegree;        
    );
    
    
    monomialDegreeHeristic.targetMonomialDegree (InterpolatedComponent) := ZZ => ( ic)->
    (
        return ic#"maxDegree";        
    );
    
    monomialDegreeHeristic.targetMonomialDegree (  BlackBoxInterpolator) := ZZ => (interpolator)->
    (
        -- if we have an interpolated component at point, then return monomialDegree of that one plus 1.  
        -- nicht gut, da 
        currentMaxDegree := max (0, interpolator.maxInterpolationDegree()) ;            
        return currentMaxDegree;        
    );
    
    monomialDegreeHeristic = newClass(BasicInterpolationMonomialDegreeHeuristic, monomialDegreeHeristic);
    return monomialDegreeHeristic;
);



ConstantInterpolationMonomialDegreeHeuristic = new Type of InterpolationMonomialDegreeHeuristic;

new ConstantInterpolationMonomialDegreeHeuristic from Thing := (E,thing)->
(
    error "not supported";
);

constantInterpolationMonomialDegreeHeuristic = method();

constantInterpolationMonomialDegreeHeuristic (BlackBoxInterpolator, ZZ) := 
    ConstantInterpolationMonomialDegreeHeuristic=>(interpolatorP, monomialDegree)->
(
    interpolator := interpolatorP;
    monomialDegreeHeristic := new MutableHashTable;
    
    monomialDegreeHeristic.targetMonomialDegree = method();
     
    monomialDegreeHeristic.targetMonomialDegree ( Matrix ) := ZZ => 
        ( point )->
    (        
         return monomialDegree;
    );
        
    monomialDegreeHeristic.targetMonomialDegree ( InterpolatedComponent ) := ZZ => 
        ( interpolatedIdeal)->
    (
        return monomialDegree;
    );
    
    monomialDegreeHeristic.targetMonomialDegree (  BlackBoxInterpolator) := ZZ => 
        ( interpolator)->
    (
         return monomialDegree;  
    );
    
    monomialDegreeHeristic = newClass(ConstantInterpolationMonomialDegreeHeuristic, monomialDegreeHeristic);
    return monomialDegreeHeristic;
)




createSimpleInterpolator = method();
createSimpleInterpolator (BlackBoxParameterSpace) := BlackBoxInterpolator => (blackBox) ->
(
    
    simpleInterpolator := MutableHashTable; -- renae to basicInterpolator?
       
    simpleInterpolator.blackBox = blackBox;
    
    monomialDegreeHeristic := null;
    jetLengthHeuristic  := null;
    
  
    
    simpleInterpolator.setJetLengthHeuristic = (jlh) ->
    (
        jetLengthHeuristic = jlh;
    );    
   
    
    simpleInterpolator.setMonomialDegreeHeristic = (mdh) ->
    (
        monomialDegreeHeristic = mdh;
    );
    
    -- I think that 'sameComponentPrecision' should coincide with interpolation jet length.
    localOnComponentPrecision := 2;

    simpleInterpolator.setOnComponentPrecision = (precision)->
    (
        localOnComponentPrecision = precision;
    );
    
    simpleInterpolator.onComponentPrecision = ()->
    (
        return localOnComponentPrecision;
    );
    
    
    -- cached component candidates
    
    -- keys are points and values are interpolated components at that points.
    --
    componentCandidatesDictionary := new MutableHashTable; -- should not be a HashTable but a type
    
    componentNamePrefix := "c";
    
    -- cached jets
    jets := new MutableHashTable;
    
  -- component id should be unique per (blackBox, interpolator)    
    nextComponentId := 0;
    
    simpleInterpolator.setComponentNamePrefix = (namePrefix)->
    (
        componentNamePrefix = namePrefix;
    );
    
    interpolationJets := new MutableHashTable;
     
    simpleInterpolator.resetInterpolation = ()->
    (
         componentCandidatesDictionary = new MutableHashTable;
         jets = new MutableHashTable; --sameComponentJets;
         interpolationJets = new MutableHashTable;
         nextComponentId = 0;
    );
    
    localAnswerStrategy := nullIfNotSmoothStrategy();
 

    simpleInterpolator.setOnComponentAnswerStrategy = method();
    simpleInterpolator.setOnComponentAnswerStrategy (Type) := Nothing => (answerStrategyType)->
    (       
        if (answerStrategyType===NullIfNotSmooth) then
        (
            localAnswerStrategy = nullIfNotSmoothStrategy();
        );
        if (answerStrategyType===ExceptionIfNotSmooth) then
        (
                localAnswerStrategy = exceptionIfNotSmooth();
        );
        if (answerStrategyType===SmoothnessInfoWithAnswerPair) then
        (
            localAnswerStrategy = smoothnessInfoWithAnswerPair();
        );      
        if (answerStrategyType===PlainTextSmoothnessInfoWithAnswerPair) then
        (
                localAnswerStrategy = plainTextSmoothnessInfoWithAnswerPair();
        );
    );
    
    simpleInterpolator.onComponentAnswerStrategies = ()->
    ( 
        return {    NullIfNotSmooth,
                    ExceptionIfNotSmooth,
                    SmoothnessInfoWithAnswerPair,
                    PlainTextSmoothnessInfoWithAnswerPair
                };
    );
    
    simpleInterpolator.setOnComponentAnswerStrategy (OnComponentAnswerStrategy) := Nothing ->(answerStrategy)->
    (
        localAnswerStrategy = answerStrategy;
    );
    
    
    simpleInterpolator.onComponentAnswerStrategy = ()->
    (
        return localAnswerStrategy  ;
    );
    
    
    -- uhh, was??? NICHT GUT? WIESO blackBox.componentCandidates()?
    -- => weil eventuell geplant war die Komponenten in der Black box zu speichern
    simpleInterpolator.componentNames = ()->
    (
        --components:= blackBox.components();
        components := values componentCandidatesDictionary;
       
        names := sort apply(components, component-> component#"name"());
        return names;        
      
    );
    
  
    simpleInterpolator.componentNameInUse = (name)->
    (
        names := simpleInterpolator.componentNames();
        return ( null=!= position(names, n->n==name) );                
    );
    
 
    
    -- question: if a name was used but is currently not, should it be forbidden to reuse it ? => no! 
    
    nextComponentName := ()->
    (        
        nextComponentId = nextComponentId+1;        
        localNextComponentName := componentNamePrefix | toString nextComponentId;        
            
        usedNames := simpleInterpolator.componentNames();
        
        while ( null=!= position(usedNames, n->n==localNextComponentName)) do
        (
            nextComponentId = nextComponentId+1;
            localNextComponentName = componentNamePrefix | toString nextComponentId;
        );      
        --print ("localNextComponentName " |localNextComponentName);
        return localNextComponentName;
    );
    
    
    
    
    simpleInterpolator.componentByName = method();
    
    simpleInterpolator.componentByName (String) := InterpolatedComponent => (name)->
    (
        selectResult := select (new HashTable from componentCandidatesDictionary, (val)->val#"name"()==name);
        if (#selectResult==0) then 
        (        
            return null;            
        );
        assert (#selectResult==1);
        return first values selectResult;        
    );        
    
    
      
    simpleInterpolator.renameComponent = method();
     
    simpleInterpolator.renameComponent (String, String) := Nothing => (name, newName)->
    (
        assert( isDerivedFrom(name, String) );
        assert( isDerivedFrom(newName, String) );
        componentToRename  := simpleInterpolator.componentByName(name);
        if (null =!= componentToRename) then 
        (
            componentToRename#"setName"(newName);
        )
        else
        (
            error ("renameComponent: no such component : '" | name  | "'" );
        );
    );
    
    simpleInterpolator.renameComponent (InterpolatedComponent, String) := Nothing => (component, newName)->
    (
        simpleInterpolator.renameComponent (component#"name"(), newName);
    );
    
    simpleInterpolator.componentNamesAt = method();
    
    simpleInterpolator.componentNamesAt ( Matrix, ZZ)  := (point, onComponentPrecisionP )->
    (
        sameComponentResult := null;
        candidates := {};
        for componentKey in keys componentCandidatesDictionary do
        (
            currentComponent := componentCandidatesDictionary#componentKey;
            sameComponentResult = catch simpleInterpolator.sameComponentAt( currentComponent, point, onComponentPrecisionP);
            if  isDerivedFrom(sameComponentResult, SingularPointException) then
            (
                return localAnswerStrategy.transformedAnswer(false,null);
            );
            if (sameComponentResult) then
            (
                  candidates = candidates | { currentComponent#"name"() };                       
            );
        );
        candidates = sort candidates;
        return localAnswerStrategy.transformedAnswer(true,candidates);        
    );
    
    
    simpleInterpolator.componentNamesAt ( Matrix)  := (point )->
    (
        return simpleInterpolator.componentNamesAt(point, localOnComponentPrecision);
    );
    
    --simpleInterpolator.componentCandidates = ()->
    --(
    --    --return new HashTable from componentCandidatesDictionary;
    --    return sort new List from values componentCandidatesDictionary;
    --);

    
    simpleInterpolator.isOnComponent = method();

    simpleInterpolator.isOnComponent ( Ideal, Matrix, ZZ) := Boolean => (componentIdeal, point, onComponentPrecisionP)->
    (
        sameComponentResult := null;
        if (sub( componentIdeal,point)!=0) then 
        (               
            return localAnswerStrategy.transformedAnswer(true, false);
        );
        sameComponentResult = catch simpleInterpolator.sameComponentAt( componentIdeal, point, onComponentPrecisionP);
        if  isDerivedFrom(sameComponentResult, SingularPointException) then
        (
            return localAnswerStrategy.transformedAnswer(false,null);            
        );  
        if  isDerivedFrom(sameComponentResult, Boolean) then
            return localAnswerStrategy.transformedAnswer(true,sameComponentResult);
        throw sameComponentResult;
    );
    
    simpleInterpolator.isOnComponent ( Ideal, Matrix) := Boolean => (componentIdeal, point)->
    (
        return simpleInterpolator.isOnComponent(componentIdeal, point);
    );
    
    simpleInterpolator.isOnComponent ( InterpolatedComponent, Matrix, ZZ) := Boolean => (componentIdeal, point, onComponentPrecisionP)->
    (
        return simpleInterpolator.isOnComponent(ideal componentIdeal, point, onComponentPrecisionP);
    );
    
    simpleInterpolator.isOnComponent ( InterpolatedComponent, Matrix) := Boolean => (componentIdeal, point)->
    (
        return simpleInterpolator.isOnComponent(ideal componentIdeal, point, localOnComponentPrecision);
    );

    
    simpleInterpolator.isOnComponent ( String, Matrix, ZZ) := Boolean => (componentName, point, onComponentPrecisionP)->
    (
       return simpleInterpolator.isOnComponent (simpleInterpolator.componentByName(componentName), point, onComponentPrecisionP);
    );
    
    simpleInterpolator.isOnComponent ( String, Matrix) := Boolean => (componentName, point)->
    (
       return simpleInterpolator.isOnComponent (simpleInterpolator.componentByName(componentName), point, localOnComponentPrecision);
    );
    
    
    simpleInterpolator.minComponentDegree = ()->
    (
        components := interpolator.components();

        if (#components>0) then 
        (
            return min apply(components, c-> max flatten degrees ideal c);
        );
        --return 0;  
        return -infinity;
    );

    simpleInterpolator.maxComponentDegree = ()->
    (
        components := interpolator.components();

        if (#components>0) then 
        (
            return max apply(components, c-> max flatten degrees ideal c);
        );
        --return 0; 
        return -infinity;
    );

    simpleInterpolator.maxInterpolationDegree = ()->
    (
        components := simpleInterpolator.components();

        if (#components>0) then 
        (
            return max apply(components, c->c#"maxDegree");
        );
        --return 0; 
        return -infinity;
    );
     
    simpleInterpolator.sameComponentAt = method();
    
    --
    -- uses a single long jet to test if a point is on a component.
    --
    
    simpleInterpolator.sameComponentAt (InterpolatedComponent, Matrix, ZZ) := Boolean  => 
        (componentIdeal, point, sameComponentPrecisionParam) ->
    (
        return simpleInterpolator.sameComponentAt (ideal componentIdeal, point, sameComponentPrecisionParam);
    );
    
    --ebenfalls cachen, da teuer. Und: mit kurzem Jet anfangen und dann immer länger werden.?
    --
    simpleInterpolator.sameComponentAt (Ideal, Matrix, ZZ) := Boolean => 
        (componentIdeal, point, sameComponentPrecisionParam)->
    (
        if (sub( componentIdeal,point)!=0) then return false;
        
        if not (jets#?point) then
        (        
            -- we have no cached jets for the given point => compute jets
            jets#point = blackBox.jetAt( point, sameComponentPrecisionParam);                        
            -- here we could also update statistic if the jet succeeded or not...
        );
       
        if length jets#point  < sameComponentPrecisionParam then 
        (
            -- we have cashed jets, but they are too short.. => compute jets of requested length
            -- improvement/optimisation: start from existing jet and enlarge it.
            jets#point = blackBox.continueJet( jets#point, sameComponentPrecisionParam);
        );
        jetP := jets#point;
        
        destEpsRng := getEpsRing(coefficientRing ring jetP#"value", sameComponentPrecisionParam);
        truncatedJetValue := sub( jetP#"value", destEpsRng );                    
        return (0 == sub( componentIdeal, truncatedJetValue ));
    );
        
    -- TODO: hide sameComponentAt !     
    simpleInterpolator.sameComponentAt (Ideal, Matrix) := Boolean => (componentIdeal, point)->
    (
        sameComponentPrecision := jetLengthHeuristic.sameComponentTargetJetLength(  componentIdeal, point);
                                                                    
        return  simpleInterpolator.sameComponentAt (componentIdeal, point,sameComponentPrecision);
    );

    interpolationJetLengthCorrection := 0 ;
    
    -- manuelle Korrektur
     simpleInterpolator.increaseInterpolationJetLength = (num)->
     (
        assert(num>0);
        interpolationJetLengthCorrection=interpolationJetLengthCorrection + num;
     );
    
    -- manuelle Korrektur
     simpleInterpolator.decreaseInterpolationJetLength = (num)->
     (
        assert(num>0);
        interpolationJetLengthCorrection=interpolationJetLengthCorrection - num;
     );
   
    -- (jk) this method should not be public. interpolate should be only provided by the black box,
    -- or even better by the interpolator.
    -- otherwise administrating component naming could get messy.
    -- so first decide: should the black box or the interpolator manage interpolated components?
    -- I vote that the interpolator should manage the components.

    interpolateAt := method();

    
   
    
    -- todo: another interpolateAt variant: pass jet(s) as parameter.
    
    cachedInterpolationJetAt := (point, targetJetLength)->
    (
      
        if not (interpolationJets#?point) then
        (        
            -- we have no cached jets for the given point => compute jets
            interpolationJets#point = blackBox.jetAt( point, targetJetLength);                        
            -- here we could also update statistic if the jet succeeded or not...
        );
       
        if length interpolationJets#point  < targetJetLength then 
        (
            -- we have cashed jets, but they are too short.. => compute jets of requested length
            -- improvement/optimisation: start from existing jet and enlarge it.
            interpolationJets#point = blackBox.continueJet( interpolationJets#point, targetJetLength);
        );        
        
        --jetP := blackBox.jetAt(point,targetJetLength);        
        jetP := interpolationJets#point;
        
        if (length interpolationJets#point  > targetJetLength) then
        (
            -- shorten the jet if necessary:
            destEpsRng := getEpsRing(coefficientRing ring jetP#"value", targetJetLength);
            truncatedJetValue := sub( jetP#"value", destEpsRng );              
            truncatedJet := jetObject(jetP#"parent", jetP#"point", truncatedJetValue, targetJetLength);
            jetP = truncatedJet;
        );
        return jetP;
    );

    interpolateAt( Matrix, ZZ, ZZ, MapHelper) := InterpolatedComponent =>
                (point,monomialMaxDegree,targetJetLength, mmap) -> 
    (
        R := mmap#"imageRing";
        mons := matrix {flatten apply(monomialMaxDegree+1, currentDegree->flatten entries basis(currentDegree,R))};
        
        -- find one jet with precision determined by jetLengthHeuristic
        
        
        -- print target Jet length at debug level?
        bblog.debug("interpolateAt targetJetLength : " |targetJetLength | " at " |toString point |" (monomialDegree="| toString monomialMaxDegree );
        
         
        
        --jetP := blackBox.jetAt(point,targetJetLength);        
        jetP := cachedInterpolationJetAt(point,targetJetLength);
        
        
        -- !!!this heuristic must be tested!!!
        -- Test: see if interpolated polynomials are in at least one
        -- irreducible component of the BlackBoxIdeal.
        jetPimage :=  (mmap#"valueAtJet")(jetP);
        --print ("jetPimage");
        --print (jetPimage);
        bareIdeal := interpolate(mons,{jetPimage});
    
        jetSet := new JetSet from jetP;
        interpolatedIdeal := createInterpolatedComponent(bareIdeal, monomialMaxDegree, jetSet, blackBox);
        return interpolatedIdeal;
    );

    interpolateAt(Matrix,ZZ,Matrix) := Ideal =>
                (point,monomialMaxDegree,mmap) -> 
    (
        targetJetLength :=  jetLengthHeuristic.interpolationTargetJetLength(  point, monomialMaxDegree);      
        targetJetLength = targetJetLength + interpolationJetLengthCorrection;
        return interpolateAt(point,monomialMaxDegree,targetJetLength, new MapHelper from mmap);
    );


    interpolateAt(Matrix,ZZ) := Ideal =>
                (point,monomialMaxDegree) -> 
    (
        targetJetLength :=  jetLengthHeuristic.interpolationTargetJetLength( point, monomialMaxDegree);      
        targetJetLength = targetJetLength + interpolationJetLengthCorrection;
        mapHelper:= createMapHelper(vars blackBox.ring, blackBox.ring);
        return interpolateAt(point,monomialMaxDegree,targetJetLength, mapHelper);
    );
    
    interpolateAt(Matrix,ZZ, ZZ) := Ideal =>
                (point,monomialMaxDegree, jetLength) -> 
    (
        mapHelper:= createMapHelper(vars blackBox.ring, blackBox.ring);
        return interpolateAt(point,monomialMaxDegree,jetLength, mapHelper);
    );
    
    
    simpleInterpolator.interpolateComponentAt = method();
    
    simpleInterpolator.interpolateComponentAt(Matrix) := Ideal => (point) -> 
    (
         localMaxMonomialDegree := max(1, simpleInterpolator.maxInterpolationDegree());
         return simpleInterpolator.interpolateComponentAt(point,localMaxMonomialDegree);
    );
    
    -- todo: rename to FORCE
    -- (JK) todo: attention: interpolateAt shares duplicate code with interpolateComponents !! 
    simpleInterpolator.interpolateComponentAt(Matrix,ZZ) := Ideal => (point, monomialMaxDegree) -> 
    (
        -- wenn schon componente da an diesem punkt, dann verbessere bis zu monomialMaxDegree
        
        --  IMPROVE AT        
        interpolationTargetDegree := null;
        interpolationJetLength := null;
        localInterpolatedIdealOrError := null;
        if (componentCandidatesDictionary#?point) then 
        (
            currentComponent := componentCandidatesDictionary#point;
            
            interpolationTargetDegree = max( monomialMaxDegree, currentComponent#"maxDegree");
           
            interpolationJetLength = jetLengthHeuristic.interpolationTargetJetLength( currentComponent, interpolationTargetDegree);            
            interpolationJetLength = interpolationJetLength + interpolationJetLengthCorrection;
            
            assert(1 == size currentComponent#"jetSet"); -- current interpolator only supports jetset with one jet
            if (currentComponent#"maxDegree" < interpolationTargetDegree or
                (length currentComponent#"jetSet"#"jets"#0) < interpolationJetLength ) then 
            (
                localInterpolatedIdealOrError = catch interpolateAt ( point,interpolationTargetDegree, interpolationJetLength );                     
            
                if (  isDerivedFrom(localInterpolatedIdealOrError,SingularPointException)) then 
                (
                    BlackBoxLogger.debug("interpolateComponents: point "| toString point| " was singular");
                    remove(componentCandidatesDictionary,point);
                )
                else  
                (   
                    -- here we do the naming magic:
                    remove(componentCandidatesDictionary,point);
                    localInterpolatedIdealOrError#"setName"( currentComponent#"name"() );
                    componentCandidatesDictionary#point = localInterpolatedIdealOrError ;  
                    return localInterpolatedIdealOrError;
                );
            )
            else
            (   
                return currentComponent;
            )                        
        )
        else
        (
            -- INTERPOLATE AT
              -- TODO: ignore incorrect points, or maybe even throw an exception or give an error
            if not blackBox.isZeroAt(point) then 
            (
                throw new PointNotOnBlackBox from {"errorMessage" => " Point " |toString point|" not on Black box"};
            );
            pointIsSingular := false;
            
            interpolationTargetDegree = monomialMaxDegree;
            -- check if component at point coincindes with one of the known components
            bSameComponentAt := false;
            for interpolatedIdeal in values componentCandidatesDictionary do
            (
                sameComponentPrecision := jetLengthHeuristic.sameComponentTargetJetLength(  interpolatedIdeal, point, interpolationTargetDegree);
                sameComponentResult := catch simpleInterpolator.sameComponentAt ( interpolatedIdeal, point, sameComponentPrecision );
                if  isDerivedFrom(sameComponentResult, SingularPointException) then
                (
                    pointIsSingular = true;
                    throw sameComponentResult;
                );
                if (sameComponentResult) then
                (
                    bSameComponentAt = true; 
                    print("warning: point is already on a component, no needs for interpolating at point");
                );                
            );
             
            -- separate compute from 
          
            interpolationJetLength = jetLengthHeuristic.interpolationTargetJetLength( point,interpolationTargetDegree);
            interpolationJetLength = interpolationJetLength + interpolationJetLengthCorrection;
            localInterpolatedIdealOrError = catch interpolateAt (point,interpolationTargetDegree, interpolationJetLength);                     
            

            if (  isDerivedFrom(localInterpolatedIdealOrError,SingularPointException)) then 
            (
                BlackBoxLogger.debug("interpolateComponents: point "| toString point| " was singular");
            )
            else  
            (                                  
                  if (  isDerivedFrom(localInterpolatedIdealOrError,InterpolatedComponent)) then
                  (
                        -- here we do the naming magic:
                        localInterpolatedIdealOrError#"setName"( nextComponentName() );
                        componentCandidatesDictionary#point = localInterpolatedIdealOrError;      
                        return localInterpolatedIdealOrError;
                  )
                  else
                  (
                     throw localInterpolatedIdealOrError;
                  );
            );
        );
    );
    
    cachedPointList := null;

    
    -- nice to have: dryRun implementieren 0> spukt alle neuen Jetlaengen 
    
    simpleInterpolator.interpolateComponentsAt = method();
    
    -- interpolateComponents should return or maintain smooth point list and eventually singular point list.

    
    simpleInterpolator.interpolateComponentsAt (List, InterpolationMonomialDegreeHeuristic, JetLengthHeuristic) := List => 
        (pointList, idh, jlh) -> 
    ( 
        if (null===pointList) then 
            error "interpolateComponents: no points given";
        
        cachedPointList = pointList;
        bblog.debug("enter InterpolateComponents");
        idealCount := 0;
        localInterpolatedIdeals := {};

        -- T := timing 
        
        err := null;
        
        -- 1. refine interpolation for known components:
        for point in keys componentCandidatesDictionary do
        (
            currentComponent := componentCandidatesDictionary#point;
            
            interpolationTargetDegree := idh.targetMonomialDegree(  currentComponent);
            interpolationJetLength := jlh.interpolationTargetJetLength( currentComponent, interpolationTargetDegree);            
            interpolationJetLength = interpolationJetLength + interpolationJetLengthCorrection;
                assert(1 == size currentComponent#"jetSet"); -- current interpolator only supports jetset with one jet
                
            if (currentComponent#"maxDegree" < interpolationTargetDegree or 
               (length currentComponent#"jetSet"#"jets"#0) < interpolationJetLength) then 
            (
                localInterpolatedIdealOrError := catch interpolateAt ( point,interpolationTargetDegree, interpolationJetLength );                     
            
                if (  isDerivedFrom(localInterpolatedIdealOrError,SingularPointException)) then 
                (
                    BlackBoxLogger.debug("interpolateComponents: point "| toString point| " was singular");
                    remove(componentCandidatesDictionary,point);
                )
                else  
                (   
                    -- here we do the naming magic:
                    remove(componentCandidatesDictionary,point);
                    localInterpolatedIdealOrError#"setName"( currentComponent#"name"() );
                    componentCandidatesDictionary#point = localInterpolatedIdealOrError ;  
                );
            )
        );
        
        for point in pointList do
        (
            BlackBoxLogger.debug("interpolateComponents: point "| toString point );

            -- TODO: ignore incorrect points, or maybe even throw an exception or give an error
            if not blackBox.isZeroAt(point) then 
            (
               continue;
            );
            pointIsSingular := false;

            -- check if component at point coincindes with one of the known components
            bSameComponentAt := false;
            for interpolatedIdeal in values componentCandidatesDictionary do
            (
                sameComponentPrecision := jlh.sameComponentTargetJetLength(  interpolatedIdeal, point);
                sameComponentResult := catch simpleInterpolator.sameComponentAt ( interpolatedIdeal, point, sameComponentPrecision );
                if  isDerivedFrom(sameComponentResult, SingularPointException) then
                (
                    pointIsSingular = true;
                    break;
                );
                if (sameComponentResult) then
                (
                    bSameComponentAt = true; 
                    --print "bSameComponentAt";
                    break;
                );                
            );
            if (bSameComponentAt) then 
            (
                BlackBoxLogger.debug("interpolateComponents: point "| toString point| " is on component!");
                continue;
            );
            if ( pointIsSingular) then 
            (
                BlackBoxLogger.debug("interpolateComponents: point "| toString point| " was singular");
                continue;
            );

            -- separate compute from 
            interpolationTargetDegree := idh.targetMonomialDegree(  point);
            interpolationJetLength := jlh.interpolationTargetJetLength( point,interpolationTargetDegree);
            interpolationJetLength = interpolationJetLength + interpolationJetLengthCorrection;
            localInterpolatedIdealOrError := catch interpolateAt (point,interpolationTargetDegree, interpolationJetLength);                     
            

            if (  isDerivedFrom(localInterpolatedIdealOrError,SingularPointException)) then 
            (
                BlackBoxLogger.debug("interpolateComponents: point "| toString point| " was singular");
            )
            else  
            (                                  
                  if (  isDerivedFrom(localInterpolatedIdealOrError,InterpolatedComponent)) then
                  (
                        -- here we do the naming magic:
                        localInterpolatedIdealOrError#"setName"( nextComponentName() );
                        componentCandidatesDictionary#point = localInterpolatedIdealOrError;                        
                  )
                  else
                  (
                     throw localInterpolatedIdealOrError;
                  );
            );
        );
        -- print "timing for loop", T#0;
                      
        return simpleInterpolator.components();
    );
    
    -- now, what should happen here?
    simpleInterpolator.interpolateComponentsAt  (List) := List=> (pointList) -> 
    (       
        -- TODO wenn kein Monomgrad angegeben ist,fange bei max (1, maxInterpolationDegree()) an
        -- und höre auf, sobald alle (glatten) Punkte auf Komponenten liegen:     
        --
        --error ("not implemented yet");        
        interpolationMaxDegree := max(1, simpleInterpolator.maxInterpolationDegree());
        return simpleInterpolator.interpolateComponentsAt(pointList,                                                        
                                                        constantInterpolationMonomialDegreeHeuristic(simpleInterpolator,
                                                                                                     interpolationMaxDegree),
                                                        jetLengthHeuristic);
    );
    
    simpleInterpolator.refineInterpolation = () -> 
    (           
        interpolationMaxDegree := 1+max(1, simpleInterpolator.maxInterpolationDegree());
        return simpleInterpolator.interpolateComponentsAt(cachedPointList,                                                        
                                                        constantInterpolationMonomialDegreeHeuristic(simpleInterpolator,
                                                                                                     interpolationMaxDegree),
                                                        jetLengthHeuristic);
    );
    
    
    simpleInterpolator.interpolateComponentsAt  (List, ZZ) := List=> (pointList, interpolationMaxDegree) -> 
    (       
        return simpleInterpolator.interpolateComponentsAt(pointList,   
                                                        constantInterpolationMonomialDegreeHeuristic(simpleInterpolator,
                                                                                                     interpolationMaxDegree),
                                                        jetLengthHeuristic);
    );
    
    
    simpleInterpolator.interpolateComponentsAt  ( ZZ) := List=> ( interpolationMaxDegree) -> 
    (
        return simpleInterpolator.interpolateComponentsAt(cachedPointList, 
                                                        constantInterpolationMonomialDegreeHeuristic(simpleInterpolator,
                                                                                                     interpolationMaxDegree),
                                                        jetLengthHeuristic
                                                        );
    );
    

    --simpleInterpolator.refineInterpolation  = () -> 
    -- (
    --    currentMaxDegree := monomialDegreeHeristic.targetMonomialDegree( blackBox, point, simpleInterpolator);
    --    return simpleInterpolator.interpolateComponents(cachedPointList, interpolationMaxDegree, sameComponentPrecision);
    --);
    
    -- simpleInterpolator.refineInterpolation  = (pointList) -> 
    --(
    --    currentMaxDegree := monomialDegreeHeristic.targetMonomialDegree( blackBox, point, simpleInterpolator);
    --    return simpleInterpolator.interpolateComponents(pointList, interpolationMaxDegree, sameComponentPrecision);
    --); 
    
    -- several components  because the point could be singular
    
    
     simpleInterpolator.componentNamesAt = method();
    
    simpleInterpolator.componentNamesAt ( Matrix, ZZ)  := (point, onComponentPrecisionP )->
    (
        sameComponentResult := null;
        candidates := {};
        for componentKey in keys componentCandidatesDictionary do
        (
            currentComponent := componentCandidatesDictionary#componentKey;
            sameComponentResult = catch simpleInterpolator.sameComponentAt( currentComponent, point, onComponentPrecisionP);
            if  isDerivedFrom(sameComponentResult, SingularPointException) then
            (
                return localAnswerStrategy.transformedAnswer(false,null);                       
            );
            if (sameComponentResult) then
            (
                  candidates = candidates | { currentComponent#"name"() };                       
            );
        );
        candidates = sort candidates;        
        return localAnswerStrategy.transformedAnswer(true,candidates);                       
    );
    
    
    simpleInterpolator.componentNamesAt ( Matrix)  := (point )->
    (
        return simpleInterpolator.componentNamesAt(point, localOnComponentPrecision);
    );
    
    
    simpleInterpolator.componentsAt = method();
    simpleInterpolator.componentsAt (Matrix, ZZ) := Thing =>( point, onComponentPrecisionP)->
    (
        sameComponentResult := null;
        candidates := {};
        for componentKey in keys componentCandidatesDictionary do
        (
            currentComponent := componentCandidatesDictionary#componentKey;
            sameComponentResult = catch simpleInterpolator.sameComponentAt( currentComponent, point, onComponentPrecisionP);
            if  isDerivedFrom(sameComponentResult, SingularPointException) then
            (
                return localAnswerStrategy.transformedAnswer(false,null);                       
            );
            if ( sameComponentResult ) then
            (
                candidates = candidates | { componentCandidatesDictionary#componentKey };
            );
        );
        candidates = sort candidates;
        return localAnswerStrategy.transformedAnswer(true,candidates);       
    );
    
    simpleInterpolator.componentsAt (Matrix) := Thing =>( point)->
    (
        return simpleInterpolator.componentsAt(point, localOnComponentPrecision);
    );
    
    simpleInterpolator.components = ( )->
    (       
        return sort values componentCandidatesDictionary;
    );
    
    simpleInterpolator.dropComponent = method();
    
    simpleInterpolator.dropComponent (InterpolatedComponent) := Nothing => (component)->
    (
        selectResult := select (new HashTable from componentCandidatesDictionary, (val)->val===component);
        if (#selectResult>0) then 
        (
            remove(componentCandidatesDictionary, first keys selectResult);
        )
        else
        (
            error ("no such component");
        );
        return;
    );
    
    
    simpleInterpolator.dropComponent (String) := Nothing => (componentName)->
    (
        selectResult := select (new HashTable from componentCandidatesDictionary, (val)->val#"name"()===componentName);
        if (#selectResult>0) then 
        (
            remove(componentCandidatesDictionary, first keys selectResult);
        )
        else
        (
            error ("no such component");
        );
        
        return;
    );
    
    
    simpleInterpolator.clearCache = ()->
    (
         jets = new MutableHashTable;
         cachedPointList = null;
    );
    
    
    simpleInterpolator = newClass(BlackBoxInterpolator, simpleInterpolator);
    
   
    -- setMonomialDegreeHeristic geht erst hier, da newClass das Objekt simpleInterpolator verändert!
    -- setJetLengthHeuristic geht erst hier, da newClass das Objekt simpleInterpolator verändert!
    simpleInterpolator.setMonomialDegreeHeristic( basicInterpolationMonomialDegreeHeuristic(simpleInterpolator));
    simpleInterpolator.setJetLengthHeuristic( basicJetLengthHeuristic(simpleInterpolator));
    return simpleInterpolator;
);

TEST /// 
    -- bug: a component name is used again after renaming.
    -- fix: always increase nextComponentId by one. 
    -- loadPackage "BlackBoxIdeals"
    errorDepth=2
    K = ZZ/11;
    R = K[x,y]
    I = ideal (x-y^2)*x;
    bbI = new BlackBoxIdeal from I;
    p1 = matrix {{0,1_K}}
    p2 = matrix {{1,1_K}}
    c1 = bbI.interpolateComponentAt(p1)
    c1Name = c1#"name"()
    bbI.renameInterpolatedComponent(c1Name, "renamedComponent")
    c2 = bbI.interpolateComponentAt(p2)
    c2Name =  c2#"name"()
    assert(c1Name != c2Name)
    
///

TEST ///
    -- bug: jet length too short 
    -- loadPackage "BlackBoxIdeals"
    kk = ZZ
    R = kk[x,y]
    I  = ideal (x*y*(x^2-y^2)*(y^4-3*x-7))
    bb = new BlackBoxIdeal from I; 
    
    p3 = matrix {{3,2_kk}}
    origin = matrix {{0, 0_kk}}
    singularPoint  = matrix{{0,0_kk}}
    bb.valuesAt p3
    jlh = basicJetLengthHeuristic(bb.interpolator);    
    maxDegree = 4;  
    mons = matrix {flatten apply(maxDegree+1, currentDegree->flatten entries basis(currentDegree,bb.ring))};
    
    jetLength = jlh.interpolationTargetJetLength(maxDegree)
    assert (jetLength >= rank source mons +10);
    
    pointList = {p1,singularPoint,p2,p3}
    pointList = {p1,p2,p3}
    pointList = {p3}
    maxDegree=1
    forcedInterpolationPrecision = 11
    iiList1 =  bb.interpolateComponentsAt(pointList,maxDegree)
    apply(iiList1, ic-> assert(1 == ic#"maxDegree"))
    maxDegree = 2
    iiList2 = bb.interpolateComponentsAt(pointList,maxDegree)   
    maxDegree=1
    iiList2b = bb.interpolateComponentsAt(pointList,1)
    apply(iiList2b, ic-> assert(2 == ic#"maxDegree"))
    
    bb.resetInterpolation()
    maxDegree = 2
    iiList2b = bb.interpolateComponentsAt(pointList,maxDegree)
    -- test zu lang !!
    --maxDegree = 5
    --bb.interpolateComponentsAt(pointList,maxDegree)
    --maxDegree = 6
    --bb.interpolateComponentsAt(pointList,maxDegree) 
    --bb.resetInterpolation()
    --maxDegree = 4
    --bb.interpolateComponentsAt(pointList,maxDegree) --ok
    --maxDegree = 5
///


TEST ///
   -- test for bug where onComponentPrecision is not respected (a longer jet was used instead of a short one)
    -- restart
    -- loadPackage "BlackBoxIdeals"
 
    setRandomSeed(42); -- ensure reproducibility in test
    K = ZZ/7;
    R = K[x,y,z,w] ;     
    line = ideal (x,y);
    conic = ideal (w,x^2+y^2-z^2);
    bbI = blackBoxIdeal intersect(line,conic);
    assert (bbI.onComponentPrecision()==2); --default onComponent precision
    pointOnLine = matrix{{0,0,1,2_K}};
    pointOnConic = matrix{{3,4,5,0_K}};
    bbI.isZeroAt(pointOnLine);
    bbI.isZeroAt(pointOnConic);
    bbI.interpolateComponentAt(pointOnLine,1);
    bbI.renameInterpolatedComponent("c1","line");
    bbI.interpolateComponentAt(pointOnConic,1);
    bbI.renameInterpolatedComponent("c2","conic");
    pointOnLineAndPlane = matrix{{0,0,1,0_K}};
    assert(bbI.isZeroAt(pointOnLineAndPlane));
    componentNames = bbI.interpolatedComponentNamesAt(pointOnLineAndPlane)    ;
    assert(1 == #componentNames);
    bbI.setOnComponentPrecision(0);
    -- now we should get both components:
    componentNames = bbI.interpolatedComponentNamesAt(pointOnLineAndPlane);
    assert(2 == #componentNames);
    assert (null =!= position(componentNames, name->name == "line"));
    assert (null =!= position(componentNames, name->name == "conic"));
    

///


TEST ///
    -- bug: jet length too short (incorrect jet length heuristic)
    -- loadPackage "BlackBoxIdeals"
    kk = QQ
    R = kk[x,y]
    I  = ideal (x*y*(x^2-y^2)*(y^4-3*x-7))
    bb = new BlackBoxIdeal from I;
    p1 = matrix {{1,1_kk}}
    p2 = matrix {{1,0_kk}}
    p3 = matrix {{3,2_kk}}
    origin = matrix {{0, 0_kk}}
    singularPoint  = matrix{{0,0_kk}}
    pointList = {p1,singularPoint,p2,p3}
    maxDegree = 1
    iiList1 = bb.interpolateComponentsAt(pointList,maxDegree)
    
    
    -- test zu lang!!!
    --maxDegree = 3
    --iiList3 = bb.interpolateComponentsAt(pointList,maxDegree) 
    --assert(3 == #iiList3);
    --bb.resetInterpolation()
    --maxDegree = 4
    --iiList4 = bb.interpolateComponentsAt(pointList,maxDegree)
    --assert(3 == #iiList4)
    
    --c3 = bb.componentsAt(p3)
    --assert (#c3 ==1)
    
    --c3 = first c3
    --assert (1==#(flatten entries gens ideal  c3))
    --assert (4==first flatten degrees ideal  c3)

    -- dies geht natuerlich nur wenn man alle Punkte angibt)
    --resultIdeal = product( apply( iiList4, ic->ideal ic))  
    -- assert (gens radical resultIdeal%  radical I
  
///


new BlackBoxInterpolator  from BlackBoxParameterSpace := (E, blackBox )->
(
    return createSimpleInterpolator(blackBox);
);
 
 
 

doc ///
   Key
        interpolate
   Headline
        find polynomials containing a list of jets
   Usage   
        I = interpolate(mons,jetList)
   Inputs  
        mons:Matrix 
            of monomials
        jetList:List
            of jets    
   Description
        Text
            Finds those linear combinations of monomials that vanish
            on the given list of jets.
               
            Lets consider a black box that describes
            a line and a plane intersecting at the origin:    
        Example      
            K = ZZ/5
            R = K[x,y,z]
            I = ideal (x*z,y*z)
            bb = blackBoxIdeal I;       
        Text
            \break 
            Consider a point on the line:
        Example
            point = matrix{{0,0,1_K}}
        Text
            \break
            and a jet of lenght 3 starting at this point and
            lying on the variety described by the black box ideal
        Example
            j = bb.jetAt(point,4)     
        Text
            \break
            Now find linear polynomials containing this jet:
        Example
            interpolate(matrix{{x,y,z}},{j})   
        Text
            Notice that polynomials containig the line are found.
            The surface is invisible to the interpolation.   
   Caveat
       This function does not estimate the lenght of the jet needed to
       get a useful answer. (The jet should be at least as long as the
       number of monomials). This is done by @TO interpolateAt @. 
   SeeAlso
       interpolateAt
///


 
doc ///
    Key
        "interpolateAt"           
    Headline
        find polynomials containing a list of jets
    Usage   
        I = interpolateAt(point,maxDegree)
        I = interpolateAt(point,maxDegree,map)
    Inputs  
        maxDegree:ZZ 
            the maximal degree of polynomials considered
        BlackBox:BlackBoxIdeal
        point: Matrix
            a point where the Blackbox vanishes    
        map: MapHelper
            
    Description
        Text
           Finds all polynomials of degree at most maxDegree
           that contain the component on which the point lies.
           If the point is not smooth, an error will be produced.
       
           Lets consider a black box that describes
           a line and a plane intersecting at the origin:
        Example      
           K = ZZ/5
           R = K[x,y,z]
           I = ideal (x*z,y*z)
           bb = blackBoxIdeal I;       
        Text
           \break 
           Consider two points on the variety described 
           by the blackbox:
        Example
           pointOnLine = matrix{{0,0,1_K}}
           pointOnPlane = matrix{{0,1,0_K}}
        Text
           \break
           Now find linear equations containing the respective
           components on which the points lie:
        Text
           The following Example has a problem:
           bb.interpolateAt(pointOnLine, 1, 10)
           bb.interpolateAt(pointOnPlane, 1, 10)
        Text
           \break
           Finding points on the different components can be done
           by running an  Experiment. Interpolating a component
           for all points found can be done by ....   
    Caveat
        This function does not work with multigraded rings.
        At the moment this has to be done by hand with @TO interpolate @. 
      
        At the moment the interpolation is done by producing one
        jet of the appropriate length. Often one could interpolate
        much faster if several shorter jets were used. (Most of the
        time is used when producing the jets)
    SeeAlso
        interpolate
        createInterpolatedComponent    
///
