@call "%~dp0PyMOLWin.exe" +2 %*
@if %ERRORLEVEL% NEQ 0 (
  @echo PyMOLWin terminated with error code %ERRORLEVEL%
  @cmd /K
)
