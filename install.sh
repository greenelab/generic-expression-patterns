#!/bin/bash
conda env create --file environment_linux.yml
conda activate generic_expression_linux
pip install -e .