
clc
clear all
close all

addpath('maxflow_mwrapper', 'maxflow_mwrapper/maxflow-v3.0');

seg_mask = segObject2D_InteractiveGraphcuts();
