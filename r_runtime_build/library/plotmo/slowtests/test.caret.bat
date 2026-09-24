@rem test.caret.bat: test plotmo on caret models
@rem Stephen Milborrow, Shrewsbury Aug 2016

@echo test.caret.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.caret.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.caret.Rout:
@echo.
@tail test.caret.Rout
@echo test.caret.R
@exit /B 1
:good1
mks.diff test.caret.Rout test.caret.Rout.save
@if %errorlevel% equ 0 goto good2:
@echo === Files are different ===
@diffps -s Rplots.ps ..\..\.#\test-reference\test.caret.save.ps
@exit /B 1
:good2
@rem test.caret.save.ps is too big to be included in the release
@rem so it is stored elsewhere
diffps Rplots.ps ..\..\.#\test-reference\test.caret.save.ps
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.caret.Rout
@rm -f Rplots.ps
@exit /B 0
