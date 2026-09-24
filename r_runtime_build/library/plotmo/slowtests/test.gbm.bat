@rem test.gbm.bat: gbm tests for plotmo and plotres

@echo test.gbm.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.gbm.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.gbm.Rout:
@echo.
@tail test.gbm.Rout
@echo test.gbm.R
@exit /B 1
:good1
mks.diff test.gbm.Rout test.gbm.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.gbm.save.ps
@exit /B 1
:good2
@rem test.gbm.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.gbm.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.gbm.Rout
@rm -f Rplots.ps
@exit /B 0
