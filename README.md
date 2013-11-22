### FiniteFieldExperiments.M2

A Macaulay2-Framework for finite field experiments for explicit and implicitly given ideals and parameter spaces
( FiniteFieldExperiments.m2)


The package provides also lifting of (polynomial) system solutions over a prime field to an extension field of rationals
( padicLift.m2 )


To use the packages,
install them first:


- start M2 (available at http://www.math.uiuc.edu/Macaulay2/Downloads/)
- run
 
```
uninstallPackage"M2Logging"
installPackage"M2Logging"
check M2Logging

uninstallPackage"IntervalPkg"
installPackage"IntervalPkg"
check IntervalPkg

uninstallPackage"BlackBoxIdeals"
installPackage"BlackBoxIdeals"
check BlackBoxIdeals

uninstallPackage"FiniteFieldExperiments"
installPackage"FiniteFieldExperiments"
check FiniteFieldExperiments

uninstallPackage"padicLift"
installPackage"padicLift"
check padicLift

```


Then look at the examples in the *experiments*-folder
or for help . e.g. by 
`viewHelp FiniteFieldExperiments`


