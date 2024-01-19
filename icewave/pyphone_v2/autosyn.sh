#!/usr/bin/env bash
time (echo $EPOCHREALTIME && adb -s 192.168.112.100:5500 shell /system/bin/time -p echo \$EPOCHREALTIME && echo $EPOCHREALTIME)>>timelog_2/t0.txt 2>&1
time (echo $EPOCHREALTIME && adb -s 192.168.112.102:5502 shell /system/bin/time -p echo \$EPOCHREALTIME && echo $EPOCHREALTIME)>>timelog_2/t2.txt 2>&1
