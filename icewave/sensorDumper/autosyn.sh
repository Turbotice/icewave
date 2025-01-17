#!/usr/bin/env bash
(gdate +%s.%N && time (echo $EPOCHREALTIME && adb -s 192.168.214.133 shell /system/bin/time -p echo \$EPOCHREALTIME && echo $EPOCHREALTIME) && gdate +%s.%N)>timelog/t0.txt 2>&1