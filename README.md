# Simulation Study Code for "Monitoring Color Image Processes"

This repository contains the code for the simulation study presented in the paper titled "Monitoring Color Image Processes" by Yarema Okhrin, Viktoriia Petruk, and Wolfgang Schmid.

## Paper Description

In this paper, we consider the problem of monitoring color image processes with spatially correlated pixels. The objective is to detect a change in the color or intensity of some pixels as soon as possible after its occurrence. We explore both the direct monitoring of color images using RGB or HSI channels and the monitoring of the corresponding grey-scale images. For dimension reduction purposes, we aggregate the pixel information over regions of interest while taking the spatial dependence into account. Three proposed GLR-type control charts are compared in an extensive simulation study.

## Important Notes

- **Not an R Package:** This code is not an R package and may require modifications for specific use cases. It is shared to make the process of the simulation study transparent.
- **Code Attribution:** Part of the code used in this simulation study is based on work by Okhrin Y, Schmid W, and Semeniuk I, from their paper "New approaches for monitoring image data" (IEEE Transactions on Image Processing 30:921–933, 2020).

## Getting Started

When examining the files in this project, begin with `calculations before sim.R`. This script is essential for setting up the simulations and includes:

- Opening a nominal image used in the study.
- Introducing functions to transform the image into a matrix of mean values over regions of interest (ROIs), to creating a covariance matrix for simulations, transform the image into gray-scale etc.
- Needed calculations (using introduced functions)
- Creating different changes for the simulation.

## Simulation Phase

During the simulation phase, you can modify the code in files related to Average Run Lengths (ARL) calculation to:

- Estimate Control Limits without introducing a change.
- Apply different change scenarios to compare the response of the control statistics.

## Authors

- Yarema Okhrin
- Viktoriia Petruk
- Wolfgang Schmid

## License

