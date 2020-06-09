@echo upgrade pymol to the latest available version
@set PREFIX=%~dp0.

"%PREFIX%\Scripts\conda" install ^
  -p "%PREFIX%" ^
  -c schrodinger pymol

@pause
