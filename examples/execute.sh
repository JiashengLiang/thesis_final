#!/bin/bash
#delect the log from last simulation 
rm iapws_log.txt ideal_log.txt
dmd isentropic_nozzle.d IAPWS_local.d 
./isentropic_nozzle