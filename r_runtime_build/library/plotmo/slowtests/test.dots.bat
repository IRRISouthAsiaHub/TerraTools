@rem test.dots.R:  test handling of dots arguments

@echo test.dots.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.dots.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.dots.Rout:
@echo.
@tail test.dots.Rout
@echo test.dots.R
@exit /B 1
:good1
mks.diff test.dots.Rout test.dots.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.dots.save.ps
@exit /B 1
:good2
@rem test.dots.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.dots.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.dots.Rout
@rm -f Rplots.ps
@exit /B 0
