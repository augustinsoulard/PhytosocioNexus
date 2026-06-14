@echo off
setlocal enabledelayedexpansion

:: Vérifie que Git LFS est installé
where git-lfs >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo ❌ Erreur : Git LFS n'est pas installé. Installez-le via https://git-lfs.com/
    pause
    exit /b 1
)

:: Vérifie qu'on est dans un dépôt Git
git rev-parse --is-inside-work-tree >nul 2>&1
if %ERRORLEVEL% neq 0 (
    echo ❌ Erreur : Ce n'est pas un dépôt Git. Exécutez ce script depuis la racine de votre dépôt.
    pause
    exit /b 1
)

:: Vérifie s'il y a des PDF modifiés
git status --porcelain | findstr /C:".pdf" >nul
if %ERRORLEVEL% equ 0 (
    echo 🔍 Nouveaux PDF détectés :
    git status --porcelain | findstr /C:".pdf"

    :: Ajoute les PDF à Git LFS
    echo ✅ Ajout des PDF à Git LFS...
    git lfs track "*.pdf"

    :: Ajoute tous les fichiers
    git add .gitattributes
    git add .

    :: Demande un message de commit
    set /p commit_message="📝 Entrez un message de commit (appuyez sur Entrée pour utiliser 'Ajout de nouveaux PDF') : "
    if "%commit_message%"=="" set commit_message=Ajout de nouveaux PDF

    :: Commit et push
    echo 🚀 Commit et push...
    git commit -m "%commit_message%"
    git push origin HEAD

    echo ✅ Terminé ! Tous les PDF ont été ajoutés et poussés.
) else (
    echo ℹ️ Aucun nouveau PDF détecté. Rien à faire.
)

pause