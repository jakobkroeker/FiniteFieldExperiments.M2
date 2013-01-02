
-- todo: root computing interface  in extra package file ? 
-- todo: parametrize path to binary
-- todo: check if gp is present/install/compile
-- optional todo: incorporate alternative root finding algorithms in GAP
callGP = method();
callGP (String,Boolean) := String=> (command,bQuiet) -> (
    progname := gppath | "gp ";
    if (bQuiet) then 
    progname = progname | " -q"; 
      gp := openInOut concatenate("!",progname," 2>/dev/null");
    --openOutAppend MEM << "cIC2 "<<get befehl<<close;
        gp  << command;
     gp  << closeOut;
     result:=get gp;
    return result;
);

-- call Pari's gp with a gp-readable command
callGP (String):= String=> (command)->
(
    quiet := true;
    return callGP( command,quiet );
)

-- given an univariate polynomial and decimal precision, creates a gp command to compute its roots with requested precision.
prepareGPRootFindingCommand = (polynomial,decimals)->
(
   command:= "/* pari/gp input: */\n";
    command = command |"\\p " | decimals |";\n";
    --command =command |"\n polynomial=" |toExternalString(polynomial) |";\n";   
    polynomialString := toExternalString(polynomial);
    polynomialString1 :=replace("ii","I",polynomialString);
    polynomialString2 :=replace( "e","E",polynomialString1);
    polynomialString3 :=replace( "p" | toString (precision polynomial) ,"",polynomialString2);
    command =command |"\n polynomial=" | polynomialString3 |";\n";
    command =command |"computedRoots=polroots(polynomial);\n";
    command =command |" modifiedroots=vector(length(computedRoots)*2);\n";
    command =command |" for(X=1,length(computedRoots),modifiedroots[(X-1)*2+1]=real(computedRoots[X]) );\n";
    command =command |" for(X=1,length(computedRoots),modifiedroots[(X-1)*2+2]=imag(computedRoots[X]) );\n";
    command=command |"print (\"beginOutput\");";
    command =command |" print (modifiedroots );";
    command =  command | "\n/* end gp input: */";
    return command;
)


-- extractArrayString; internal. Extracts from a string the bracketed part. 
-- used to convert pari's root computation output to Macaulay.
-- precondition : '[' and ']' occur only once. 
--
-- Example:  extractArrayString("bla[1,s],bka") returns "[1,s]"
extractArrayString = method()
extractArrayString (String) :=String=>(str)->
(
    -- todo: count '[' and ']' should be 1 .
     pos1 := regex("\\[",str);
     pos2 := regex("\\]",str);
     if pos1===null or pos2===null then error "could not extract the array...";
     pos1 = pos1#0#0+1;
     pos2 = pos2#0#0;
     assert( pos2>pos1 );
     subs:=substring(pos1,pos2-pos1,str);
    return subs;
)

testExtractArrayString =()->
(
    assert(extractArrayString("bla[1,s],bka")=="1,s");
    --replace("I","ii",str)
)

-- convert an pari root computation result (for the prepareGPRootFindingCommand)  from the output  string to Macaulay2
convertPariComplexNumArrayToMacaulay = ( pariArrayStr, decimalPrecision )->
(
    -- count '[' and ']' should be 1 .
    --subs :=     extractArrayString(pariArrayStr );

    localstr  := replace("I","ii",pariArrayStr);
    localstr1 := replace("E","e",localstr);
    localstrList := separate(",",localstr1);
    localstrListStr := "";
 
    apply(localstrList,el-> (if (localstrListStr!="")then  localstrListStr=localstrListStr |","; localstrListStr=localstrListStr | el ));
    --
    bitprecision := round( decimalPrecision*log(2,10) );
    --
   --error in /tmp/tmp.iJnwE3
   -- Macaulay2-usage from GAP regularily caused an error. writing the input to a file "stringRootArray.m2" was an (failed) attempt to debug the problem.
   --ofile := openOut ("stringRootArray.m2");
   --ofile << "localstrListStr :=" << (toExternalString localstrListStr) << ";\n";
   --ofile << "decimalPrecision :=" << decimalPrecision << ";\n";
   --ofile << "bitprecision := " << bitprecision << ";\n";
   --ofile << "-- toExternalString (round( decimalPrecision*log(2,10) )) << \n";
   --ofile << "bitprecision :=" <<  toExternalString (round( decimalPrecision*log(2,10) )) << ";\n";
   --ofile << "log2:= " << (log(2,10)) << ";\n";
   localstrListStr = for el in separate(",",localstrListStr) list
    (
          el1:=replace(" ","",el);
          el2:=replace("\n","",el1);
          el3:=replace("\t","",el2);
          pos1:=regex("e",el3);
          el4:=null;
          if (pos1=!=null) then 
            el4=replace("e",  ("p"|bitprecision |"e"),el3)
          else
            el4=el3 | "p"|bitprecision;
         el4
    );
    ----complexRootArray := apply(((#localstrListStr)//2),pos->(print (localstrListStr#(pos*2));print  (localstrListStr#(pos*2+1) );) );
    ----print localstrListStr;
    --todo: aus einem toExternalString p.. entfernen! - warum?  todo still 
    --ofile << (toExternalString localstrListStr) << close;
    complexRootArray := apply(((#localstrListStr)//2),pos->((value (localstrListStr#(pos*2)))+ii* (value (localstrListStr#(pos*2+1)))));
    return complexRootArray;
)


computeRootsWithGP=method(Options=>{"verbose"=>false});
computeRootsWithGP (RingElement,ZZ):=opts->(polynomial,decimalPrecision)->
(
    try (char ring polynomial) then (assert (char ring polynomial==0)) else () ;
    rootFindingCommand := prepareGPRootFindingCommand(polynomial,decimalPrecision);
    --print rootFindingCommand;
    resultStr := callGP(rootFindingCommand);
    if (opts#"verbose") then 
    (
        print "resultStr";
        print resultStr;
    );
    pos1 := regex("beginOutput",resultStr);
    assert(pos1=!=null);
    subs := substring(pos1#0#0,resultStr);
    --
    arrayStr := extractArrayString(subs);
    rootList := convertPariComplexNumArrayToMacaulay(arrayStr,decimalPrecision);
    return rootList;
)

doc ///
    Key
        computeRootsWithGP        
        (computeRootsWithGP, RingElement, ZZ)
    
    Headline
        compute univariate polynomials complex roots using Pari's implementation
    Usage   
        computeRootsWithGP( univariatePolynomial, decimalPrecision)
    Inputs  
        univariatePolynomial:RingElement
             univariate polynomial
        decimalPrecision:ZZ
             decimal precision 
    Outputs
        : List
             a complex root list for the input polynomial
    Description
        Example          
        Text
           \break  Compute roots for  (x^2 + 1)
        Example          
            rng := CC[x];
    	    polynomial := x^2+1;   decimalPrecision := 12;
            rootList := computeRootsWithGP(polynomial,decimalPrecision)
        Text
           \break  Compute roots for  (x^2+1/3)
        Example          
            rng = QQ[x];
    	    polynomial =  x^2+1/3; 
            rootList = computeRootsWithGP(polynomial,decimalPrecision)
///

-- testing the (ugly) interface to pari' root computation
testComputeRootsWithGP=()->
(

    -- testing  root computation for a polynomial with complex coefficients
    x := null;
    rng := CC[(symbol x)];
    x = (gens rng)#0;
    polynomial := x^2+1;
    decimalPrecision := 10;
    rootList := computeRootsWithGP(polynomial,decimalPrecision);
    assert(#rootList==2);
    apply(rootList, el-> assert(sub(polynomial,x=>el)==0));

    -- testing  root computation for apolynomial with rational coefficients
    x = null;
    rng = QQ[(symbol x)];
    x = (gens rng)#0;
    polynomial = x^2+1/3;
    decimalPrecision = 10;
    rootList = computeRootsWithGP(polynomial,decimalPrecision);
    assert(#rootList==2);
    apply(rootList, el-> assert(sub(polynomial,x=>el)==0));

    -- testing root computation for a polynomial with indexed variable.
    y := null;  y   = symbol y;
    rng = QQ[(y_1)];
    yy := (gens rng)#0;
    polynomial = yy^2+1/3;
    decimalPrecision = 10;
    rootList = computeRootsWithGP(polynomial,decimalPrecision);
    assert(#rootList==2);
    apply(rootList, el-> assert(sub(polynomial,y_1=>el)==0));

)

-- test, if Macaulay2 is able to call Pari's gp binary.
pariGpIsPresent=()->
(
    try (
    bQuiet:=false;
    result := callGP("",bQuiet);
    pos1:=regex("GP/PARI",result);
    return (pos1=!=null);
    )
    then ()
    else( return false;) ;
)


configurePari=()->
(
    gppath = (options padicLift).Configuration#"gppath";
    if not pariGpIsPresent() then 
        gppath = prefixDirectory | currentLayout#"programs";      
    
    --if (not pariGpIsPresent() ) then 
    --    error "please install Pari/GP or build Macaulay2 with pari (add 'pari' to the --enable-build-libraries configure parameter. ) to use this package !"        
)


configurePari();
