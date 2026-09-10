@echo off
setlocal
echo ============================================
echo  SangerScanner - Windows exe builder
echo ============================================
echo.
echo This builds in an ISOLATED virtual environment so it does not
echo pick up unrelated Anaconda packages (PyQt6 / PySide6 / etc.)
echo that break the build or bloat the exe.
echo.

REM ---- [1/4] create a clean build environment ----
echo [1/4] Creating clean build environment (build_venv)...
if exist build_venv rmdir /s /q build_venv
python -m venv build_venv
if errorlevel 1 (
    echo.
    echo Could not create venv. Make sure Python 3.9-3.12 is installed and on PATH.
    echo   Check with:  python --version
    pause
    exit /b 1
)

call build_venv\Scripts\activate.bat

REM ---- [2/4] install only what SangerScanner needs ----
echo.
echo [2/4] Installing dependencies...
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
if errorlevel 1 (
    echo.
    echo pip install failed - see the errors above.
    pause
    exit /b 1
)

REM ---- [3/4] build ----
echo.
echo [3/4] Building SangerScanner.exe with PyInstaller...
pyinstaller --onefile --windowed --name SangerScanner ^
    --hidden-import=Bio.SeqIO.AbiIO ^
    --collect-submodules=Bio ^
    --exclude-module=PyQt5 ^
    --exclude-module=PyQt6 ^
    --exclude-module=PySide2 ^
    --exclude-module=PySide6 ^
    --exclude-module=IPython ^
    --exclude-module=notebook ^
    ab1_gui.py
if errorlevel 1 (
    echo.
    echo Build failed - see the errors above.
    pause
    exit /b 1
)

REM ---- [4/4] done ----
echo.
echo ============================================
echo  Done!  dist\SangerScanner.exe
echo.
echo  That single file is the whole program - copy it
echo  anywhere and double-click to run. No Python needed.
echo ============================================
pause
endlocal
