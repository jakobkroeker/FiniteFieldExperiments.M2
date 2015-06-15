


.PHONY: installPadicLift installBlackBoxIdeals installFiniteFieldExperiments  uninstallPadicLift \
        uninstallBlackBoxIdeals uninstallFiniteFieldExperiments installM2Logging uninstallM2Logging \
        install uninstall
 
manualinstall: installPadicLift installBlackBoxIdeals installFiniteFieldExperiments installM2Logging
   

installPadicLift:  
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift.m2 ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift  ~/.Macaulay2/local/share/Macaulay2/ -R
	@echo -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/padicLift.m2 "> ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh
	@echo "rm  -rf ~/.Macaulay2/local/share/Macaulay2/padicLift/" >>  ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh


installM2Logging: 
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp M2Logging.m2 ~/.Macaulay2/local/share/Macaulay2/
	@echo  -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/M2Logging.m2 " > ~/.Macaulay2/local/share/Macaulay2/M2LoggingUninstall.sh


installBlackBoxIdeals: 
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp BlackBoxIdeals.m2 ~/.Macaulay2/local/share/Macaulay2/
	@echo  -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/BlackBoxIdeals.m2 " > ~/.Macaulay2/local/share/Macaulay2/BlackBoxIdealsUninstall.sh

installFiniteFieldExperiments: 
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp FiniteFieldExperiments.m2 ~/.Macaulay2/local/share/Macaulay2/
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperiments
	cp FiniteFieldExperiments/Interpolation.m2 ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperiments/
	@echo  -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperiments.m2 " > ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperimentsUninstall.sh
	@echo  -e "#"'!'"/bin/bash \n rm  -rf ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperiments/ " > ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperimentsUninstall.sh


uninstallPadicLift:
	bash ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh

uninstallM2Logger:
	bash ~/.Macaulay2/local/share/Macaulay2/M2LoggerUninstall.sh


uninstallBlackBoxIdeals:
	bash ~/.Macaulay2/local/share/Macaulay2/BlackBoxIdealsUninstall.sh

uninstallFiniteFieldExperiments:
	bash ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperimentsUninstall.sh


manualuninstall: uninstallPadicLift uninstallBlackBoxIdeals uninstallFiniteFieldExperiments uninstallM2Logging

check:
	M2 --script testFFexperiments.m2

#	M2 --script testPadicLift.m2

test: check

# will not work on mac (or at all? should consider auxiliary files!)
install:
	echo "path = append(path,\""`pwd`"/\")" >installPackages.m2
	cat installPackages.template >> installPackages.m2
	M2 < installPackages.m2

uninstall:
	M2 --script uninstallPackage.m2

