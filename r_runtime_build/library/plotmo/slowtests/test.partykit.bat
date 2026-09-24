@rem test.partykit.bat

@echo test.partykit.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.partykit.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.partykit.Rout:
@echo.
@tail test.partykit.Rout
@echo test.partykit.R
@exit /B 1
:good1
mks.diff test.partykit.Rout test.partykit.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.partykit.save.ps
@exit /B 1
:good2
@rem test.partykit.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.partykit.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.partykit.Rout
@rm -f Rplots.ps
@exit /B 0
