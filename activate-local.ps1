$ProjectRoot = Split-Path -Parent $MyInvocation.MyCommand.Path
$env:PYTHONPATH = "$ProjectRoot\.packages;$ProjectRoot\src"
Write-Host "DreamPep local packages enabled."
Write-Host "Python: C:\Users\86159\.cache\codex-runtimes\codex-primary-runtime\dependencies\python\python.exe"

