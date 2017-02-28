

bblog := BlackBoxLogger;
-- InterpolatedIdeal = new Type of MutableHashTable;

InterpolatedIdeal = new Type of HashTable;

new InterpolatedIdeal from Thing :=  (InterpolatedIdealAncestor,l)->
(  
    ---type check?
    error("InterpolatedIdeal from Thing not supported or not implemented");
);

InterpolatedIdeal ? InterpolatedIdeal := (i1,i2)->
(
    return (i2#"name"() ? i2#"name"());
)
--new InterpolatedIdeal from MutableHashTable :=  (InterpolatedIdealAncestor,l)->
--(  
    ----type check?
    --return l;
--);

--new InterpolatedIdeal from List :=  (InterpolatedIdealAncestor,l)->
--(  
--    return new InterpolatedIdeal from new MutableHashTable from l;
--);

net (InterpolatedIdeal) := Net =>(ii)->
(
    return netInterpolatedIdeal(ii);
)

netInterpolatedIdeal = (II)->
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

-- find polynomials containing a component
-- via interpolation

-- interpolate()
--
-- find linear combinations of monomials containing a given jet
-- it seems that currently jetList have only one entry.
--

interpolate = method();
interpolate = (mons, jetList) -> 
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


-- interpolateBB()
--
-- find all polynomials of degree smallerEqual than maxDegree containing the component of BB containing the point
-- 
-- this is the most basic simple minded implementation where only one very long jet is considered.
--

createInterpolatedIdeal = method();
createInterpolatedIdeal( Ideal, ZZ, Thing, BlackBoxParameterSpace ) := InterpolatedIdeal => 
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
    result = newClass(InterpolatedIdeal, result);
    return result;
);

createInterpolatedIdeal( Ideal, ZZ, String ) := InterpolatedIdeal => 
                       (I,maxDegree,  description)->
(
    return createInterpolatedIdeal (I,maxDegree, null, description);
);

ideal ( InterpolatedIdeal) := Ideal => (ii)->
(
    return ii#"ideal";
)

bare = method();
bare ( InterpolatedIdeal) := Ideal => (ii)->
(
    return ii#"ideal";
)

TEST ///
    -- loadPackage "BlackBoxIdeals"
    R = ZZ[x,y]
    I = ideal (x*y)
    bb = new BlackBoxIdeal from I
    p = matrix {{1,0_QQ}}
    jetset = new JetSet from jetAt(bb,p,1)
    II = createInterpolatedIdeal(I, 1, jetset, bb)
    ideal II
    bare II     
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

constantJetLengthHeuristic = (targetJetLength)->
(
    jetLengthEstimator:= new MutableHashTable;
    jetLengthEstimator.setJetLength = (newTargetJetLength)->
    (
        targetJetLength = newTargetJetLength;
    );
        
        
    jetLengthEstimator.interpolationTargetJetLength = method();
     
    -- for interpolating at a point
    --
    jetLengthEstimator.interpolationTargetJetLength ( BlackBoxParameterSpace,
                                                      BlackBoxInterpolator, 
                                                      Matrix,
                                                      ZZ                      ) :=  ZZ =>
        ( blackBox, 
          interpolator, 
          point,      
          monomialDegree )->
    (
        return targetJetLength;
    );
    
    jetLengthEstimator.interpolationTargetJetLength ( BlackBoxParameterSpace,
                                                      BlackBoxInterpolator, 
                                                      InterpolatedIdeal,
                                                      ZZ                      ) :=  ZZ =>
        ( blackBox, 
          interpolator, 
          interpolatedIdeal,      
          monomialDegree )->
    (
        return targetJetLength;
    );
    
    jetLengthEstimator.interpolationTargetJetLength ( BlackBoxParameterSpace,
                                                      BlackBoxInterpolator,                                                    
                                                      ZZ                      ) :=  ZZ =>
        ( blackBox, 
          interpolator, 
          monomialDegree )->
    (
        return targetJetLength;
    );
    
    
    jetLengthEstimator.sameComponentTargetJetLength = method();
    -- component, point (check if point on component)
    --
    jetLengthEstimator.sameComponentTargetJetLength (BlackBoxParameterSpace, 
                                                     BlackBoxInterpolator, 
                                                     InterpolatedIdeal, 
                                                     Matrix                )     :=  ZZ => 
        (blackBox, 
         interpolator, 
         component, 
         point)->
    (
        return targetJetLength;
    );
    
    result := newClass(ConstantJetLengthHeuristic, jetLengthEstimator);
    return result;
);



BasicJetLengthHeuristic = new Type of JetLengthHeuristic;

new BasicJetLengthHeuristic from Thing := (E,thing)->
(
    error "not supported";
);

basicJetLengthHeuristic = ()->
(
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
    
     jetLengthEstimator.interpolationTargetJetLength (BlackBoxParameterSpace,BlackBoxInterpolator,ZZ) :=  ZZ =>
        (blackBox, interpolator, monomialDegree)->
    (
        return rank super basis(monomialDegree, blackBox.coefficientRing) + localAdditionalJetLength;
    );
    
     
    jetLengthEstimator.interpolationTargetJetLength (BlackBoxParameterSpace,BlackBoxInterpolator, Matrix,ZZ) :=  ZZ =>
        (blackBox, interpolator, point, monomialDegree)->
    (
        return rank super basis(monomialDegree, blackBox.coefficientRing) + localAdditionalJetLength;
    );
    
    jetLengthEstimator.sameComponentTargetJetLength = method();
    
     -- component, point (check if point on component)
    jetLengthEstimator.sameComponentTargetJetLength (BlackBoxParameterSpace, BlackBoxInterpolator, InterpolatedIdeal, Matrix) :=  ZZ => 
        (blackBox, interpolator, component, point)->
    (
        monomialDegree := component#"maxDegree";        
        return  rank super basis(monomialDegree, blackBox.coefficientRing) + localAdditionalJetLength;
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

basicInterpolationMonomialDegreeHeuristic = ()->
(
    monomialDegreeHeristic := new MutableHashTable;
    
    monomialDegreeHeristic.targetMonomialDegree = method();
    
    
    monomialDegreeHeristic.targetMonomialDegree (BlackBoxParameterSpace, BlackBoxInterpolator, Matrix ) := ZZ => 
        (blackBox,interpolator, point )->
    (        
        currentMaxDegree := 0;
        
        components := interpolator.componentsAt(point);
        if (#components>0) then 
        (
               currentMaxDegree = max (0, apply (components, c-> max degrees flatten ideal c)) ;  
        );               
        return currentMaxDegree+1;        
    );
    
    
    monomialDegreeHeristic.targetMonomialDegree (BlackBoxParameterSpace, BlackBoxInterpolator,InterpolatedIdeal) := ZZ => 
        (blackBox, interpolator, ic)->
    (
        return ic#"maxDegree"+1;        
    );
    
    monomialDegreeHeristic.targetMonomialDegree (BlackBoxParameterSpace,  BlackBoxInterpolator) := ZZ => 
        (blackBox,  interpolator)->
    (
        -- if we have an interpolated component at point, then return monomialDegree of that one plus 1.  
        -- nicht gut, da 
        currentMaxDegree := max (0, interpolator.maxInterpolationDegree()) ;            
        return currentMaxDegree+1;        
    );
    
    monomialDegreeHeristic = newClass(BasicInterpolationMonomialDegreeHeuristic, monomialDegreeHeristic);
    return monomialDegreeHeristic;
)

ConstantInterpolationMonomialDegreeHeuristic = new Type of InterpolationMonomialDegreeHeuristic;

new ConstantInterpolationMonomialDegreeHeuristic from Thing := (E,thing)->
(
    error "not supported";
);

constantInterpolationMonomialDegreeHeuristic = (monomialDegree)->
(
    monomialDegreeHeristic := new MutableHashTable;
    
    monomialDegreeHeristic.targetMonomialDegree = method();
     
    monomialDegreeHeristic.targetMonomialDegree (BlackBoxParameterSpace, BlackBoxInterpolator, Matrix ) := ZZ => 
        (blackBox,interpolator, point )->
    (        
         return monomialDegree;
    );
        
    monomialDegreeHeristic.targetMonomialDegree (BlackBoxParameterSpace,BlackBoxInterpolator, InterpolatedIdeal ) := ZZ => 
        (blackBox, ic, interpolator)->
    (
        return monomialDegree;
    );
    
    monomialDegreeHeristic.targetMonomialDegree (BlackBoxParameterSpace,  BlackBoxInterpolator) := ZZ => 
        (blackBox,  interpolator)->
    (
         return monomialDegree;  
    );
    
    monomialDegreeHeristic = newClass(ConstantInterpolationMonomialDegreeHeuristic, monomialDegreeHeristic);
    return monomialDegreeHeristic;
)

BlackBoxInterpolator = new Type of HashTable;


createSimpleInterpolator = method();
createSimpleInterpolator (BlackBoxParameterSpace) := BlackBoxInterpolator => (blackBox) ->
(
    
    simpleInterpolator := MutableHashTable; -- renae to basicInterpolator?
       
    simpleInterpolator.blackBox = blackBox;
    
    jetLengthHeuristic := basicJetLengthHeuristic();
    
    
    simpleInterpolator.setJetLengthHeuristic = (jlh) ->
    (
        jetLengthHeuristic = jlh;
    );    
    
    
    monomialDegreeHeristic := basicInterpolationMonomialDegreeHeuristic();
    
    
    simpleInterpolator.setMonomialDegreeHeristic = (mdh) ->
    (
        monomialDegreeHeristic = mdh;
    );
    
    -- I think that 'sameComponentPrecision' should coincide with interpolation jet length.
    --sameComponentPrecision := 5;

    --simpleInterpolator.setSameComponentPrecision = (precision)->
    --(
    --    sameComponentPrecision=precision;
    --);
    
    
    -- cached component candidates
    
    -- keys are points and values are interpolated components at that points.
    --
    componentCandidatesDictionary := new MutableHashTable; -- should not be a HashTable but a type
    
    componentNamePrefix := "c";
    
    
    simpleInterpolator.setComponentNamePrefix = (namePrefix)->
    (
        componentNamePrefix = namePrefix;
    );
    
    -- component id should be unique per (blackBox, interpolator)    
    nextComponentId := 1;
    
    
    -- uhh, was??? NICHT GUT? WIESO blackBox.componentCandidates()?
    -- => weil eventuell geplant war die Komponenten in der Black box zu speichern
    simpleInterpolator.componentNames = ()->
    (
        --components:= blackBox.componentCandidates();
        components := values componentCandidatesDictionary;
       
        names := apply(components, component-> component#"name"());
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
    
    
    -- (jk) well, I think we do not need setComponentName(point, newName) , 
    -- and, this could be unfortunate if the point is not smooth
    -- and we have several components there
    --
    
    simpleInterpolator.componentByName = method();
    
    simpleInterpolator.componentByName (String) := InterpolatedIdeal => (name)->
    (
        selectResult := select (new HashTable from componentCandidatesDictionary, (val)->val#"name"()==name);
        if (#selectResult==0) then 
        (        
            return null;            
        );
        assert (#selectResult==1);
        return first values selectResult;        
    );        
    
    
      
    
    simpleInterpolator.renameComponent = (name, newName)->
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
    
    
    simpleInterpolator.componentNamesAt = (point)->
    (
        candidates := {};
        for componentKey in keys componentCandidatesDictionary do
        (
            currentComponent := componentCandidatesDictionary#componentKey;
            if ( simpleInterpolator.isOnComponent( currentComponent, point)) then
            (
                candidates = candidates | { currentComponent#"name"() };
            );
        );
        return   candidates;        
    );
    
    -- cached jets

    jets := new MutableHashTable;

    simpleInterpolator.componentCandidates = ()->
    (
        --return new HashTable from componentCandidatesDictionary;
        return sort new List from values componentCandidatesDictionary;
    );

  

    -- todo: rename to isProbablyOnComponent
    --
    simpleInterpolator.isOnComponent = method();

    simpleInterpolator.isOnComponent ( Ideal, Matrix) := Boolean => (componentIdeal, point)->
    (
        return (sub( componentIdeal,point)==0);
    );
    
    simpleInterpolator.isOnComponent ( InterpolatedIdeal, Matrix) := Boolean => (componentIdeal, point)->
    (
        return simpleInterpolator.isOnComponent(ideal componentIdeal, point);
    );
    
    simpleInterpolator.isOnComponent ( String, Matrix) := Boolean => (componentName, point)->
    (
       return simpleInterpolator.isOnComponent (simpleInterpolator.componentByName(componentName), point);
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
        components := interpolator.components();

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
    
    simpleInterpolator.sameComponentAt (InterpolatedIdeal, Matrix, ZZ) := Boolean  => 
        (componentIdeal, point, onComponentPrecisionParam) ->
    (
        return simpleInterpolator.sameComponentAt (ideal componentIdeal, point, onComponentPrecisionParam);
    );
    
    simpleInterpolator.sameComponentAt (Ideal, Matrix, ZZ) := Boolean => 
        (componentIdeal, point, sameComponentPrecisionParam)->
    (
        if (sub( componentIdeal,point)!=0) then return false;
        
        bblog.debug ("enter simpleInterpolator.sameComponentAt");
        if not (jets#?point) then
        (        
            -- we have no cached jets for the given point => compute jets
            jets#point = jetAt( blackBox ,point, sameComponentPrecisionParam);
            -- here we could also update statistic if the jet succeeded or not...
        );
        jetP := jets#point;  
        
        if (jetP === null) then
        (
            bblog.debug ("leave 1 simpleInterpolator.sameComponentAt");
            errorMessage:= {"errorMessage"=>"jets not succeeded. Point "|  toString point | "is not smooth?";};
            throw new SingularPointException from errorMessage
        );

        if jetP#"jetLength" < sameComponentPrecisionParam then 
        (
            -- we have cashed jets, but they are too short.. => compute jets of requested length
            -- improvement/optimisation: start from existing jet and enlarge it.
            jets#point = continueJet( blackBox ,jetP, sameComponentPrecisionParam);
            jetP = jets#point;  
        );

        --if jetP#"succeeded" then 
        if (jetP =!= null) then
        (
            -- check if the point is on the component or not
            -- print "succeeded";            
            -- testing throw new SingularPointException from {"errorMessage"=>"jets not succeeded. Point is not smooth?";};
            bblog.debug ("leave 3 simpleInterpolator.sameComponentAt");
            return (0 == sub( componentIdeal, jetP#"value" ));
        )
        else
        (
            bblog.debug ("leave 2 simpleInterpolator.sameComponentAt");
            throw new SingularPointException from {"errorMessage"=>"jets not succeeded. Point is not smooth?";};
        );

    );
    
    simpleInterpolator.sameComponentAt (Ideal, Matrix) := Boolean => (componentIdeal, point)->
    (
        sameComponentPrecision := jetLengthHeuristic.sameComponentTargetJetLength(blackBox,
                                                                            simpleInterpolator,
                                                                            componentIdeal,
                                                                            point);
                                                                    
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

    simpleInterpolator.interpolateAt = method();

    -- todo: another interpolateAt variant: pass jet(s) as parameter.

    simpleInterpolator.interpolateAt( Matrix, ZZ, ZZ, MapHelper) := InterpolatedIdeal =>
                (point,monomialMaxDegree,targetJetLength, mmap) -> 
    (
        R := mmap#"imageRing";
        mons := matrix {flatten apply(monomialMaxDegree+1, currentDegree->flatten entries basis(currentDegree,R))};
        
        -- find one jet with precision determined by jetLengthHeuristic
        
        
        -- print target Jet length at debug level?
        bblog.debug("interpolateAt targetJetLength : " |targetJetLength | " at " |toString point |" (monomialDegree="| toString monomialMaxDegree );
        
        jetP := blackBox.jetAt(point,targetJetLength);        
        -- !!!this heuristic must be tested!!!
        -- Test: see if interpolated polynomials are in at least one
        -- irreducible component of the BlackBoxIdeal.
        jetPimage :=  (mmap#"valueAtJet")(jetP);
        bareIdeal:= interpolate(mons,{jetPimage});
    
        jetSet := new JetSet from jetP;
        interpolatedIdeal := createInterpolatedIdeal(bareIdeal, monomialMaxDegree, jetSet, blackBox);
        return interpolatedIdeal;
    );

    simpleInterpolator.interpolateAt(Matrix,ZZ,Matrix) := Ideal =>
                (point,monomialMaxDegree,mmap) -> 
    (
        targetJetLength :=  jetLengthHeuristic.interpolationTargetJetLength(blackBox, simpleInterpolator, point, monomialMaxDegree);      
        targetJetLength = targetJetLength + interpolationJetLengthCorrection;
        return simpleInterpolator.interpolateAt(point,monomialMaxDegree,targetJetLength, new MapHelper from mmap);
    );


    simpleInterpolator.interpolateAt(Matrix,ZZ) := Ideal =>
                (point,monomialMaxDegree) -> 
    (
        targetJetLength :=  jetLengthHeuristic.interpolationTargetJetLength(blackBox, simpleInterpolator, point, monomialMaxDegree);      
        targetJetLength = targetJetLength + interpolationJetLengthCorrection;
        return simpleInterpolator.interpolateAt(point,monomialMaxDegree,targetJetLength, createMapHelper(vars blackBox.ring, blackBox.ring));
    );
    
    simpleInterpolator.interpolateAt(Matrix,ZZ, ZZ) := Ideal =>
                (point,monomialMaxDegree, jetLength) -> 
    (
        return simpleInterpolator.interpolateAt(point,monomialMaxDegree,jetLength,createMapHelper(vars blackBox.ring, blackBox.ring));
    );
    
    cachedPointList := null;

    
    -- nice to have: dryRun implementieren 0> spukt alle neuen Jetlaengen 
    
    simpleInterpolator.interpolateComponents = method();
    
    -- interpolateComponents should return or maintain smooth point list and eventually singular point list.
    
    simpleInterpolator.interpolateComponents (List, ZZ, ZZ) := List => 
        (pointList, interpolationMaxDegree, targetJetLength) -> 
    (
        -- hier bauen wir konstante Heuristiken fr maxMonomialDegree und sameComponentPrecisionParam:
        jlh := constantJetLengthHeuristic(targetJetLength);
        
        mdh := constantInterpolationMonomialDegreeHeuristic(interpolationMaxDegree);
        return simpleInterpolator.interpolateComponents(pointList,mdh, jlh);
    );
    
    
    
    simpleInterpolator.interpolateComponents (List, InterpolationMonomialDegreeHeuristic, JetLengthHeuristic) := List => 
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
            
            interpolationTargetDegree := idh.targetMonomialDegree(blackBox,simpleInterpolator, currentComponent);
            interpolationJetLength := jlh.interpolationTargetJetLength(blackBox,simpleInterpolator,currentComponent, interpolationTargetDegree);
            interpolationJetLength = interpolationJetLength + interpolationJetLengthCorrection;
            if (currentComponent#"maxDegree" < interpolationTargetDegree) then 
            (
                localInterpolatedIdealOrError := catch simpleInterpolator.interpolateAt ( point,interpolationTargetDegree, interpolationJetLength );                     
            
                if (  isDerivedFrom(localInterpolatedIdealOrError,SingularPointException)) then 
                (
                    BlackBoxLogger.debug("interpolateComponents: point "| toString point| " was singular");
                    remove(componentCandidatesDictionary,point);
                )
                else  
                (   
                    -- here we do the naming magic:
                    localInterpolatedIdealOrError#"setName"( nextComponentName() );
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
                sameComponentPrecision := jlh.sameComponentTargetJetLength(blackBox,simpleInterpolator, interpolatedIdeal, point);
                sameComponentResult := simpleInterpolator.sameComponentAt ( interpolatedIdeal, point, sameComponentPrecision );
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
            interpolationTargetDegree := idh.targetMonomialDegree(blackBox,simpleInterpolator, point);
            interpolationJetLength := jlh.interpolationTargetJetLength(blackBox,simpleInterpolator,point,interpolationTargetDegree);
            interpolationJetLength = interpolationJetLength + interpolationJetLengthCorrection;
            localInterpolatedIdealOrError := catch simpleInterpolator.interpolateAt (point,interpolationTargetDegree, interpolationJetLength);                     
            

            if (  isDerivedFrom(localInterpolatedIdealOrError,SingularPointException)) then 
            (
                BlackBoxLogger.debug("interpolateComponents: point "| toString point| " was singular");
            )
            else  
            (                                  
                  if (  isDerivedFrom(localInterpolatedIdealOrError,InterpolatedIdeal)) then
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
                      
        return simpleInterpolator.componentCandidates();
    );
     simpleInterpolator.interpolateComponents  (List) := List=> (pointList) -> 
    (       
        return simpleInterpolator.interpolateComponents(pointList,   
                                                        monomialDegreeHeristic,
                                                        jetLengthHeuristic);
    );
    simpleInterpolator.interpolateComponents  (List, ZZ) := List=> (pointList, interpolationMaxDegree) -> 
    (       
        return simpleInterpolator.interpolateComponents(pointList,   
                                                        constantInterpolationMonomialDegreeHeuristic(interpolationMaxDegree),
                                                        jetLengthHeuristic);
    );
    
    simpleInterpolator.interpolateComponents  ( ZZ) := List=> ( interpolationMaxDegree) -> 
    (
        return simpleInterpolator.interpolateComponents(cachedPointList, 
                                                        constantInterpolationMonomialDegreeHeuristic(interpolationMaxDegree),
                                                        jetLengthHeuristic
                                                        );
    );
    
    simpleInterpolator.interpolateComponents  ( ZZ,ZZ) := List=> ( interpolationMaxDegree, sameComponentPrecisionParam) -> 
    (
        return simpleInterpolator.interpolateComponents(cachedPointList, 
                                                        constantInterpolationMonomialDegreeHeuristic(interpolationMaxDegree),
                                                        constantJetLengthHeuristic(sameComponentPrecisionParam)
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
    simpleInterpolator.componentsAt = ( point)->
    (
        candidates := {};
        for componentKey in keys componentCandidatesDictionary do
        (
            if ( simpleInterpolator.isOnComponent( componentCandidatesDictionary#componentKey, point)) then
            (
                candidates = candidates | { componentCandidatesDictionary#componentKey };
            );
        );
        return candidates;
    );
    
    simpleInterpolator.dropComponent= method();
    
    simpleInterpolator.dropComponent (InterpolatedIdeal) := Nothing => (component)->
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
     
    return simpleInterpolator;
);



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
            j = jetAt(bb,point,4)     
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
       createInterpolatedIdeal
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
        Example
           bb.interpolateAt(pointOnLine, 1,10)
           bb.interpolateAt(pointOnPlane, 1,10)
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
        createInterpolatedIdeal    
///
