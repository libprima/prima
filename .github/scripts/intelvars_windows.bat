for /f "tokens=* usebackq" %%f in (`dir /b "C:\Program Files (x86)\Intel\oneAPI\compiler\" ^| findstr /V latest ^| sort`) do @set "LATEST_VERSION=%%f"
@call "C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\env\vars.bat"

echo %LATEST_VERSION%
where ifort.exe
where ifx.exe

echo ONEAPI_VER=%LATEST_VERSION%>> %GITHUB_ENV%
echo ONEAPI_ROOT=C:\Program Files (x86)\Intel\oneAPI>> %GITHUB_ENV%
echo IFORT_COMPILER21=C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\windows>> %GITHUB_ENV%
echo PATH=%PATH%;C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\windows\bin\intel64;C:\Program Files (x86)\Intel\oneAPI\compiler\%LATEST_VERSION%\windows\bin>> %GITHUB_ENV%
