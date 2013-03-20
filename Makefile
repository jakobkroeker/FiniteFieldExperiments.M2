
.PHONY: installPadicLift installIdealBlackBoxes uninstallPadicLift uninstallIdealBlackBoxes  install uninstall
 
install: installPadicLift installIdealBlackBoxes
   

installPadicLift: 
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift.m2 ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift  ~/.Macaulay2/local/share/Macaulay2/ -R
	@echo -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/padicLift.m2 "> ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh
	@echo "rm  -rf ~/.Macaulay2/local/share/Macaulay2/padicLift/" >>  ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh

installIdealBlackBoxes: 
	mkdir -p   ~/.Macaulay2/local/share/Macaulay2/
	cp idealBlackBoxes.m2 ~/.Macaulay2/local/share/Macaulay2/
	cp padicLift  ~/.Macaulay2/local/share/Macaulay2/ -R
	@echo  -e "#"'!'"/bin/bash \n rm  ~/.Macaulay2/local/share/Macaulay2/idealBlackBoxes.m2 " > ~/.Macaulay2/local/share/Macaulay2/idealBlackBoxesUninstall.sh


uninstallPadicLift:
	bash ~/.Macaulay2/local/share/Macaulay2/padicLiftUninstall.sh


uninstallIdealBlackBoxes:
	bash ~/.Macaulay2/local/share/Macaulay2/idealBlackBoxesUninstall.sh


uninstall: uninstallPadicLift uninstallIdealBlackBoxes
