@rem test.modguide.bat: test model1 and model2 (linmod examples) in modguide.pdf

@echo test.modguide.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.modguide.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.modguide.Rout:
@echo.
@tail test.modguide.Rout
@echo test.modguide.R
@exit /B 1
:good1
mks.diff test.modguide.Rout test.modguide.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.modguide.save.ps
@exit /B 1
:good2
@rem test.modguide.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.modguide.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.modguide.Rout
@rm -f Rplots.ps
@exit /B 0
