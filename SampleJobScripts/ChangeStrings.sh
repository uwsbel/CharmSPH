#!/bin/bash

## Change Strings

sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes1/N1/*

sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes2/N1/*
sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes2/N2/*
sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes2/N4/*

sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes4/N1/*
sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes4/N2/*
sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes4/N4/*

sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes8/N1/*
sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes8/N2/*
sed -i 's/H=0.025/H=0.05/g' Euler_Runs_h0.05/numNodes8/N4/*

