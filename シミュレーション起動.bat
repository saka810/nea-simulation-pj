@echo off
rem ===================================================================
rem  幾何音響シミュレーション  起動用（ダブルクリックで開く）
rem
rem  ・この .bat はリポジトリ直下に置いたまま使う（場所を基準に動く）
rem  ・プロジェクトフォルダを、このファイルの上にドラッグ＆ドロップすると
rem    そのプロジェクトを開いた状態で始まる
rem  ・黒い画面（実行ログ）は閉じないこと。計算の進み具合が出る
rem ===================================================================
setlocal
title 幾何音響シミュレーション

set "HERE=%~dp0"
set "PY=%HERE%.venv\Scripts\python.exe"
set "APP=%HERE%geosim\app.py"

if not exist "%APP%" goto NOAPP
if not exist "%PY%" goto NOVENV

rem geosim はベア名インポートなので、そこを作業フォルダにして実行する
cd /d "%HERE%geosim"
"%PY%" "%APP%" %*
if errorlevel 1 goto FAILED

endlocal
exit /b 0


:NOAPP
echo.
echo  [エラー] geosim\app.py が見つかりません。
echo.
echo   この起動ファイルは、リポジトリ（nea-simulation-pj）の直下に
echo   置いたまま使ってください。別の場所へコピーすると動きません。
echo.
echo   いまの場所: %HERE%
echo.
pause
endlocal
exit /b 1


:NOVENV
echo.
echo  [エラー] 実行環境（.venv）が見つかりません。
echo.
echo   .venv は端末ごとに作り直すものです（Git では配られません）。
echo   PowerShell を開いて、下の 3 行を順に実行してください。
echo.
echo     cd "%HERE%"
echo     py -3.10 -m venv .venv
echo     .venv\Scripts\python -m pip install -r requirements.txt
echo.
echo   ※ Python は 3.10 を明示すること（チーム方針。.python-version 参照）
echo.
pause
endlocal
exit /b 1


:FAILED
echo.
echo  --------------------------------------------------------------
echo   プログラムが途中で終了しました。
echo   上に出ているメッセージが原因です（そのまま報告してください）。
echo  --------------------------------------------------------------
echo.
pause
endlocal
exit /b 1
