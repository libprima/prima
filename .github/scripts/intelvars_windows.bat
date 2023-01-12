for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

echo %LATEST_VERSION%
where ifort.exe
where ifx.exe

echo PATH=%PATH%;C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\bin\intel64;C:\Program Files (x86)\Intel\oneAPI\compiler\latest\windows\bin >> %GITHUB_ENV%
