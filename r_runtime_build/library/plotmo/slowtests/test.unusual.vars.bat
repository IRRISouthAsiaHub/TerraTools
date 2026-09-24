@rem test.unusual.vars.bat: test unusual variable names and formulas
@rem                        this file was first created for plotmo 3.6.0 (Sep 2020)

@echo test.unusual.vars.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.unusual.vars.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.unusual.vars.Rout:
@echo.
@tail test.unusual.vars.Rout
@echo test.unusual.vars.R
@exit /B 1
:good1
mks.diff test.unusual.vars.Rout test.unusual.vars.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.unusual.vars.save.ps
@exit /B 1
:good2
@rem test.unusual.vars.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.unusual.vars.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.unusual.vars.Rout
@rm -f Rplots.ps
@exit /B 0
