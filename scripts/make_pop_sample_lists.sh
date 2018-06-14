#!/bin/sh

# ----------------------------------------------------------------------------------------
# --- Make population sample lists
# ----------------------------------------------------------------------------------------

grep "eRHG" data/East_AGR_RHG.txt | cut -f1 > data/eRHG_IDs.txt
grep "eAGR" data/East_AGR_RHG.txt | cut -f1 > data/eAGR_IDs.txt
