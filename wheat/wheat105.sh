#!/bin/sh


runseries ()
{
  lsysBasename="${1}"
  for d in ${dayList} ; do
    ltransgl -geometry 800x600+0+0 -d $d -p 0,20,40 -r 0,90,0 -i ${lsysBasename}_d${d}.ppm ${lsysBasename}.trl 
  done
}

dayList='000 001 010 020 026 027 032 033 034 035 036 037 038 040 041 042 043 044 045 046 047 048 049 050 051 052 053 054 055 056 057 058 059 060 061 062 063 064 068 069 070 071 072 073 074 075 076 077 078 079 080 081 082 083 084 085 086 087 088 089 092 093 094 095 096 097 099 100 101 102 103 104 105 106 108 109 110'
runseries wheat105

# for wheatModel in singleshoot alwaysgreen senescence ; do
#  runseries wheat105_${wheatModel}
# done

