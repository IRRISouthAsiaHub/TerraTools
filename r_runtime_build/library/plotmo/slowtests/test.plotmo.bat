@rem test.plotmo.bat: this does a regression test of plotmo
@rem Stephen Milborrow Apr 2007 Petaluma

@echo test.plotmo.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.plotmo.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.plotmo.Rout:
@echo.
@tail test.plotmo.Rout
@echo test.plotmo.R
@exit /B 1
:good1
mks.diff test.plotmo.Rout test.plotmo.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.plotmo.save.ps
@exit /B 1
:good2
@rem test.plotmo.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.plotmo.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.plotmo.Rout
@rm -f Rplots.ps
@exit /B 0
