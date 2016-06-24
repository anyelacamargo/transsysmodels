#!/bin/sh


runseries ()
{
  lsysBasename="${1}"
  for d in ${dayList} ; do
    ltransgl -geometry 800x600+0+0 -d $d -p 0,20,40 -r 0,90,0 -i ${lsysBasename}_d${d}.ppm ${lsysBasename}.trl 
  done
}

dayList='000 001 010 020 030 040 050 060 070 080 090 100 110'
runseries wheat105

# for wheatModel in singleshoot alwaysgreen senescence ; do
#  runseries wheat105_${wheatModel}
# done

