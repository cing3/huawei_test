@echo off
setlocal enabledelayedexpansion

set SOLUTION_FILE=%1
set INPUT_FILE=%2
set OUTPUT_FILE=%3

if "%SOLUTION_FILE%"=="" (
    echo Usage: %0 ^<solution_file^> ^<input_file^> ^<output_file^>
    echo Example: %0 .\demos\Solution.cpp .\data\sample.in .\data\sample.out
    echo Supported languages: C, CPP, JAVA, PYTHON
    exit /b 5
)
if "%INPUT_FILE%"=="" (
    echo Usage: %0 ^<solution_file^> ^<input_file^> ^<output_file^>
    echo Example: %0 .\demos\Solution.cpp .\data\sample.in .\data\sample.out
    echo Supported languages: C, CPP, JAVA, PYTHON
    exit /b 5
)
if "%OUTPUT_FILE%"=="" (
    echo Usage: %0 ^<solution_file^> ^<input_file^> ^<output_file^>
    echo Example: %0 .\demos\Solution.cpp .\data\sample.in .\data\sample.out
    echo Supported languages: C, CPP, JAVA, PYTHON
    exit /b 5
)

set "LANG_TYPE="

for %%A in ("%SOLUTION_FILE%") do set "EXT=%%~xA"

if /i "%EXT%"==".c" set LANG_TYPE=C
if /i "%EXT%"==".cpp" set LANG_TYPE=CPP
if /i "%EXT%"==".java" set LANG_TYPE=JAVA
if /i "%EXT%"==".py" set LANG_TYPE=PYTHON

if "%LANG_TYPE%"=="" (
    echo Unsupported solution file extension: %SOLUTION_FILE%
    exit /b 5
)

for %%f in ("%SOLUTION_FILE%") do (
    set "BASE_NAME=%%~nf"
    set "SOLUTION_PATH=%%~dpf"
)

if /i "%LANG_TYPE%"=="C" goto compile_cpp
if /i "%LANG_TYPE%"=="CPP" goto compile_cpp
if /i "%LANG_TYPE%"=="JAVA" goto compile_java
if /i "%LANG_TYPE%"=="PYTHON" goto compile_python

echo Unsupported language: %LANG_TYPE%
exit /b 5

:compile_cpp
echo Compiling Solution.cpp...
g++ -std=c++11 %SOLUTION_FILE% -o "%SOLUTION_PATH%\%BASE_NAME%.exe" -O2 -lpthread
if errorlevel 1 (
    echo Compilation error: %SOLUTION_FILE% compilation failed
    exit /b 2
)
if not exist "%SOLUTION_PATH%\%BASE_NAME%.exe" (
    echo File not found: "%SOLUTION_PATH%\%BASE_NAME%.exe" file does not exist
    exit /b 6
)
goto run_runner

:compile_java
echo Compiling Solution.java...
javac -d %SOLUTION_PATH% %SOLUTION_FILE%
if errorlevel 1 (
    echo Compilation error: %SOLUTION_FILE% compilation failed
    exit /b 2
)
if not exist "%SOLUTION_PATH%\%BASE_NAME%.class" (
    echo File not found: "%SOLUTION_PATH%\%BASE_NAME%.class" file does not exist
    exit /b 6
)
goto run_runner

:compile_python
echo Python does not need to be compiled
goto run_runner

:run_runner
echo Begin running ...
Runner.exe %INPUT_FILE% %OUTPUT_FILE% %LANG_TYPE% %SOLUTION_PATH% %BASE_NAME%

echo Finish process
endlocal
exit /b 0
