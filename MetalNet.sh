#!/bin/bash
input_ID=1535637622
python predict.py $input_ID
python network_filter.py $input_ID
