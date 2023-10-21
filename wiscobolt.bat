@echo off
for /f %%a in ('powershell -command "Get-Date -Format MM-dd-yy_HH-mm-ss"') do set datetime=%%a

set folder_name=Results\%datetime%
set "folder_name=%folder_name:\=/%"  REM Replace backslashes with forward slashes

if not exist "%folder_name%" mkdir "%folder_name%"

@echo on
make all

@echo off
set temp_file="%folder_name%\temp_output.txt"

wiscobolt.exe "%folder_name%" | tee %temp_file%
echo wiscobolt.exe

@echo on
make clean

@echo off
move /Y %temp_file% "%folder_name%\output.txt" > nul
