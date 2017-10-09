#!/bin/bash
#delect the log from last simulation 
rm log.txt
dmd isentropic_nozzle.d IAPWS_local.d 
./isentropic_nozzle