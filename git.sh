#!/bin/bash

pyhton RAA_eChem.Thread_MA_V_3.py

wait

git init
git add .
git commit -m "Firt commit"
git remote add origin https://github.com/Mahsa791/Test-ma.git
git push origin HEAD:master
