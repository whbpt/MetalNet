#!/bin/bash
input_ID=$1
python predict.py $input_ID
python network_filter.py $input_ID
