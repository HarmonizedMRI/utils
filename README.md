# HarmonizedMRI Utilities

Shared utility functions used across multiple repositories in the HarmonizedMRI organization.

## Scope

This repository contains reusable code that is not specific to any single MRI acquisition, reconstruction pipeline, or study. Typical examples include:

- Raw data I/O 
- Image reconstruction utilities
- General MRI processing routines

Functions in this repository are intended to be called by higher-level projects such as:

- SMS-EPI
- B0shimming
- ArbEPI 

## What does *not* belong here

This repository should **not** contain:

- Pulse sequence implementations
- Study-specific processing pipelines
- BIDS dataset construction scripts
- Scanner-specific installation workflows
- Example datasets

Those belong in the corresponding project repositories.

## Status

This repository is under active development and serves as a common code base for the HarmonizedMRI software ecosystem.

