


.PHONY: installPadicLift installBlackBoxIdeals uninstallPadicLift uninstallBlackBoxIdeals  install uninstall
 
install: installPadicLift installBlackBoxIdeals
   

installPadicLift:  
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift.m2 ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift  ~/.Macaulay2/local/share/Macaulay2/ -R
	@echo -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/padicLift.m2 "> ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh
	@echo "rm  -rf ~/.Macaulay2/local/share/Macaulay2/padicLift/" >>  ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh

installBlackBoxIdeals: 
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp BlackBoxIdeals.m2 ~/.Macaulay2/local/share/Macaulay2/
	cp FailingExport.m2 ~/.Macaulay2/local/share/Macaulay2/
	@echo  -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/BlackBoxIdeals.m2 " > ~/.Macaulay2/local/share/Macaulay2/BlackBoxIdealsUninstall.sh


uninstallPadicLift:
	bash ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh


uninstallBlackBoxIdeals:
	bash ~/.Macaulay2/local/share/Macaulay2/uninstallBlackBoxIdeals.sh


uninstall: uninstallPadicLift uninstallBlackBoxIdeals
