@rem test.plotmo3.bat: extra tests for plotmo version 3 and higher

@echo test.plotmo3.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.plotmo3.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.plotmo3.Rout:
@echo.
@tail test.plotmo3.Rout
@echo test.plotmo3.R
@exit /B 1
:good1
mks.diff test.plotmo3.Rout test.plotmo3.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.plotmo3.save.ps
@exit /B 1
:good2
@rem test.plotmo3.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.plotmo3.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.plotmo3.Rout
@rm -f Rplots.ps
@exit /B 0
