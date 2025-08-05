#!/usr/bin/env python3
"""
Check requirements for PEST Monte Carlo analysis
Run this before submitting the job to verify everything is set up correctly
"""

import os
import sys
import subprocess
from pathlib import Path

def check_python_packages():
    """Check if required Python packages are available"""
    print("Checking Python packages...")
    packages = ['pyemu', 'flopy', 'pandas', 'numpy', 'matplotlib']
    
    for package in packages:
        try:
            __import__(package)
            print(f"  ✓ {package} available")
        except ImportError:
            print(f"  ✗ {package} not available - install with: pip install {package}")

def check_executables():
    """Check if PEST++ executables are available"""
    print("\nChecking executables...")
    
    # Check in the bin directory
    bin_dir = Path("/home/tfo46/pakipaki_model01/bin")
    if not bin_dir.exists():
        print(f"  ✗ Bin directory not found: {bin_dir}")
        return
    
    executables = ['pestpp-swp.exe', 'mf6.exe', 'mp7.exe']
    
    for exe in executables:
        exe_path = bin_dir / exe
        if exe_path.exists():
            print(f"  ✓ {exe} found at {exe_path}")
            if os.access(exe_path, os.X_OK):
                print(f"    ✓ {exe} is executable")
            else:
                print(f"    ✗ {exe} is not executable - run: chmod +x {exe_path}")
        else:
            print(f"  ✗ {exe} not found in {bin_dir}")

def check_project_structure():
    """Check if the project structure matches expectations"""
    print("\nChecking project structure...")
    
    base_path = Path("/home/tfo46/pakipaki_model01")
    
    # Check critical directories
    dirs_to_check = [
        "manual_builds",
        "manual_builds/dvc_scripts",
        "bin"
    ]
    
    for dir_name in dirs_to_check:
        dir_path = base_path / dir_name
        if dir_path.exists():
            print(f"  ✓ {dir_name} directory exists")
        else:
            print(f"  ✗ {dir_name} directory missing")

def check_setup_file():
    """Check if setup.py can be imported and paths make sense"""
    print("\nChecking setup.py configuration...")
    
    try:
        # Add the script directory to Python path
        script_dir = Path("/home/tfo46/pakipaki_model01/manual_builds/dvc_scripts")
        sys.path.insert(0, str(script_dir))
        
        from setup import MODEL_NAME, BIN_DIR, TEMP_DIR, PEST_DIR
        
        print(f"  ✓ setup.py imported successfully")
        print(f"  Model name: {MODEL_NAME}")
        print(f"  Bin directory: {BIN_DIR}")
        print(f"  Template directory: {TEMP_DIR}")
        print(f"  PEST directory: {PEST_DIR}")
        
    except ImportError as e:
        print(f"  ✗ Could not import setup.py: {e}")
    except Exception as e:
        print(f"  ✗ Error with setup.py: {e}")

def main():
    print("PEST Monte Carlo Requirements Check")
    print("=" * 40)
    
    check_python_packages()
    check_executables()
    check_project_structure()
    check_setup_file()
    
    print("\n" + "=" * 40)
    print("Requirements check complete!")
    print("\nIf any items show ✗, please fix them before running the job.")
    print("To install missing Python packages, run:")
    print("  pip install pyemu flopy pandas numpy matplotlib")

if __name__ == "__main__":
    main()
