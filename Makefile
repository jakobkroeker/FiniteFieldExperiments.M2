


.PHONY: installPadicLift installBlackBoxIdeals installFiniteFieldExperiments  uninstallPadicLift uninstallBlackBoxIdeals uninstallFiniteFieldExperiments install uninstall
 
install: installPadicLift installBlackBoxIdeals installFiniteFieldExperiments
   

installPadicLift:  
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift.m2 ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift  ~/.Macaulay2/local/share/Macaulay2/ -R
	@echo -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/padicLift.m2 "> ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh
	@echo "rm  -rf ~/.Macaulay2/local/share/Macaulay2/padicLift/" >>  ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh

installBlackBoxIdeals: 
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp BlackBoxIdeals.m2 ~/.Macaulay2/local/share/Macaulay2/
	@echo  -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/BlackBoxIdeals.m2 " > ~/.Macaulay2/local/share/Macaulay2/BlackBoxIdealsUninstall.sh

installFiniteFieldExperiments: 
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp FiniteFieldExperiments.m2 ~/.Macaulay2/local/share/Macaulay2/
	@echo  -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperiments.m2 " > ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperimentsUninstall.sh


uninstallPadicLift:
	bash ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh


uninstallBlackBoxIdeals:
	bash ~/.Macaulay2/local/share/Macaulay2/BlackBoxIdealsUninstall.sh

uninstallFiniteFieldExperiments:
	bash ~/.Macaulay2/local/share/Macaulay2/FiniteFieldExperimentsUninstall.sh


uninstall: uninstallPadicLift uninstallBlackBoxIdeals uninstallFiniteFieldExperiments
