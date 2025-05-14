#!/usr/bin/env bash
sudo ifconfig eno2 192.168.0.10
sudo systemctl restart isc-dhcp-server.service
sudo systemctl status isc-dhcp-server.service
