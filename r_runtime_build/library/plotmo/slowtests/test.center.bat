@rem test.center.bat: test plotmo's center and ndiscrete args
@rem Stephen Milborrow, Berea Apr 2011

@echo test.center.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.center.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.center.Rout:
@echo.
@tail test.center.Rout
@echo test.center.R
@exit /B 1
:good1
mks.diff test.center.Rout test.center.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.center.save.ps
@exit /B 1
:good2
@rem test.center.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.center.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.center.Rout
@rm -f Rplots.ps
@exit /B 0
