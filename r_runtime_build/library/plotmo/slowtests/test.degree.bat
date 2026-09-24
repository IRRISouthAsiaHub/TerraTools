@rem test.degree.bat: test plotmo's degree1 and degree2 args with character arguments

@echo test.degree.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.degree.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.degree.Rout:
@echo.
@tail test.degree.Rout
@echo test.degree.R
@exit /B 1
:good1
mks.diff test.degree.Rout test.degree.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.degree.save.ps
@exit /B 1
:good2
@rem test.degree.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.degree.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.degree.Rout
@rm -f Rplots.ps
@exit /B 0
