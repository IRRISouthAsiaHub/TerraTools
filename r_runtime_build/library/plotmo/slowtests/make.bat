@rem plotmo/inst/slowtests/make.bat

@call test.plotmo.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.printcall.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.dots.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.plotmo.dots.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.plotmo.x.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.plotmo.args.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.degree.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.modguide.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.linmod.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.fac.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.plotmo3.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.center.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.plotres.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.partdep.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.unusual.vars.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.non.earth.bat
                                            @if %errorlevel% NEQ 0 goto err

@rem The following miscellaneous models are in alphabetical order

@call test.c50.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.caret.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.gbm.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.glmnet.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.glmnetUtils.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.mlr.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.parsnip.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.partykit.bat
                                            @if %errorlevel% NEQ 0 goto err
@call test.pre.bat
                                            @if %errorlevel% NEQ 0 goto err
@rem we also run the earth package tests in \a\r\earth\inst\slowtests\make.bat
@cd \a\r\earth\inst\slowtests
                                            @if %errorlevel% NEQ 0 goto err
@call make.bat
                                            @if %errorlevel% NEQ 0 goto err
@cd \a\r\plotmo\inst\slowtests
@goto done
:err
@echo ==== ERROR ====
:done
@exit /B  0
