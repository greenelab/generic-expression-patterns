#!/bin/bash
conda env create --file environment.yml
conda activate generic_expression
pip install -e .