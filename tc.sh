#!/bin/bash 
export LD_LIBRARY_PATH=/projappl/project_2001565/TeraChem/lib:$LD_LIBRARY_PATH
export TeraChem=/projappl/project_2001565/TeraChem 

#Cleanup Scratch directory
if [ -e scr.qm ] ; then rm -r -f scr.qm ; fi

$TeraChem/bin/terachem $1 > tc-qmmm.out
set stat=$?

cp -f scr.qm/c0 c0
cp -f scr.qm/ca0 ca0
cp -f scr.qm/cb0 cb0
cp -f scr.qm/grad.xyz grad.xyz

exit $stat

