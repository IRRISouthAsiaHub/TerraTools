@rem test.parsnip.bat
@rem Stephen Milborrow Sep 2020 Petaluma

@echo test.parsnip.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.parsnip.R
@if %errorlevel% equ 0 goto good1
@echo R returned errorlevel %errorlevel%, see test.parsnip.Rout:
@echo.
@tail test.parsnip.Rout
@echo test.parsnip.R
@exit /B 1
:good1
@rem second egrep gets rid of random messages issued by library(tidymodels)
@rem could perhaps use suppressPackageStartupMessages() instead
@egrep -v "Fit time:| Use | Dig | Learn | Search |^\* " test.parsnip.Rout >test.parsnip.Rout2
mv test.parsnip.Rout2 test.parsnip.Rout
mks.diff test.parsnip.Rout test.parsnip.Rout.save
@if %errorlevel% equ 0 goto good2
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.parsnip.save.ps
@exit /B 1
:good2
@rem test.parsnip.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.parsnip.save.ps
@if %errorlevel% equ 0 goto good3
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.parsnip.Rout test.parsnip.Rout
@rm -f Rplots.ps
@exit /B  0
