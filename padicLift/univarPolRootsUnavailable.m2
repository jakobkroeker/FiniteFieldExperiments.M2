
-- only dummy methods   -pari/gp installation is not present 

prepareGPRootFindingCommand = (polynomial,decimals)->
(
   return null;
)


extractArrayString = method()
extractArrayString (String) :=String=>(str)->
(
      return null;
)

testExtractArrayString =()->
(
      return true;
)


convertPariComplexNumArrayToMacaulay = ( pariArrayStr, decimalPrecision )->
(
       return null;
)


computeRootsWithGP=method(Options=>{"verbose"=>false});
computeRootsWithGP (RingElement,ZZ):=opts->(polynomial,decimalPrecision)->
(
         return null;
)

doc ///
    Key
        computeRootsWithGP        
        (computeRootsWithGP, RingElement, ZZ)
    Headline
        unavailable -  pari/gp required
   Description
       Text
        This method requires pari/gp which could not be found on this computer when installPackage was run.
	Install pari/gp and do
	\break
	\break uninstallPackage"padicLift"
	\break restart
	\break loadPackage"padicLift"
	\break installPackage"padicLift"

///

testComputeRootsWithGP=()->
(
    return true
)

