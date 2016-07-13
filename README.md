# Code and data for Fleming, Massoni et al. (in prep)

This repo contains Matlab analysis code and data to reproduce analyses presented in Fleming, Massoni, Gajdos & Vergnaud (in prep) "Metacognition about the past and future: Quantifying common and distinct influences on prospective and retrospective judgments of self-performance"

FlemingMassoni_NOC_analysis.m  computes various measures of prospective and retrospective accuracy, confidence (overconfidence and calibration) and metacognitive sensitivity (AUROC2 and the adjusted normalised discrimination index, ANDI). This script calls the functions brier_index and type2roc.

FlemingMassoni_Pmodel_wrapper.m fits learning models to the prospective judgments as a function of past performance or past confidence. This script calls the function fitPconf to perform gradient descent on model parameters.
