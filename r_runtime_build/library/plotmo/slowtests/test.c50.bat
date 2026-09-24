@rem test.c50.bat: c50 tests for plotmo and plotres

@echo test.c50.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.c50.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.c50.Rout:
@echo.
@tail test.c50.Rout
@echo test.c50.R
@exit /B 1
:good1
mks.diff test.c50.Rout test.c50.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.c50.save.ps
@exit /B 1
:good2
@rem test.c50.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.c50.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.c50.Rout
@rm -f Rplots.ps
@exit /B 0
