
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
