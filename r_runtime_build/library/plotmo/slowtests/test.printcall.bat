@rem test.printcall.R:  test printcall

@echo test.printcall.bat
@"C:\PROGRA~1\R\R-4.5.2\bin\x64\R.exe" CMD BATCH --quiet --vanilla test.printcall.R
@if %errorlevel% equ 0 goto good1:
@echo R returned errorlevel %errorlevel%, see test.printcall.Rout:
@echo.
@tail test.printcall.Rout
@echo test.printcall.R
@exit /B 1
:good1
mks.diff test.printcall.Rout test.printcall.Rout.save
@if %errorlevel% equ 0 goto good3:
@echo === Files are different ===
@exit /B 1
:good3
@rm -f test.printcall.Rout
@exit /B 0
