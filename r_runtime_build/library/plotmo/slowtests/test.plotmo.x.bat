@rem test.plotmo.x.bat: test plotmo_x and related functions

@echo test.plotmo.x.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.plotmo.x.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.plotmo.x.Rout:
@echo.
@tail test.plotmo.x.Rout
@echo test.plotmo.x.R
@exit /B 1
:good1
mks.diff test.plotmo.x.Rout test.plotmo.x.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.plotmo.x.save.ps
@exit /B 1
:good2
@rem test.plotmo.x.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.plotmo.x.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.plotmo.x.Rout
@rm -f Rplots.ps
@exit /B 0
