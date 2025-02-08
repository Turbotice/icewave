#!/usr/bin/env bash
adb devices > Connected
while read line; do echo $line; done < Connected
#check the length of the adb identifier, 12 : redmi10, 16:geo-phone
while read line; do id=`echo $line | awk -F ' ' '{print $1}'`; if ((${#id}==16))||((${#id}==12))||((${#id}==11)); then adress=$id; fi;done < Connected
while read line; do id=`echo $line | awk -F ' ' '{print $2}'`; if [[ $id == $adress ]]; then num=`echo $line | awk -F ' ' '{print $1}'`; fi; done < adb_usb_liste 
#if ["$id" = "$id_usb"]; then echo 'yes'; fi; done < adb_usb_liste
echo $num $adress
adb -s $adress tcpip $((5500+$num))
sleep 1
adb -s $adress connect 192.168.0.$((100+$num)):$((5500+$num))
