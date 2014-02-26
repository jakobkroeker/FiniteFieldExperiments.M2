-- testing padicLift package

uninstallPackage"padicLift"
restart
installPackage("padicLift",UserMode =>true)
check (padicLift,UserMode =>true)
loadPackage("padicLift",Reload =>true)

--viewHelp padicLift

----------
uninstallPackage"BlackBoxIdeals"
restart
loadPackage"BlackBoxIdeals"
installPackage"BlackBoxIdeals"
jetAt
