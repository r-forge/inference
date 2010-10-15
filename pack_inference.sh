#! /bin/bash

cd ~/Documents/Computing/R_packages/On_R-Forge/inference \
    && rm -rf ./pkg/man/* \
    && R CMD roxygen -d -s pkg \
    && R CMD check pkg \
    && R CMD build pkg \
    && sudo R CMD INSTALL inference*.tar.gz