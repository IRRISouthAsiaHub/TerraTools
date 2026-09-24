@rem test.linmod.bat: test example S3 model in linmod.R

@echo test.linmod.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.linmod.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.linmod.Rout:
@echo.
@tail test.linmod.Rout
@echo test.linmod.R
@exit /B 1
:good1
mks.diff test.linmod.Rout test.linmod.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.linmod.save.ps
@exit /B 1
:good2
@rem test.linmod.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.linmod.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.linmod.Rout
@rm -f Rplots.ps
@exit /B 0
