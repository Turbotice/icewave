import pprint

import requests

from cst import ADB_PORT, BASH_SCRIPTS_PATH, IP_BASE, PHYPHOX_PORT, DUMP, ID_ADB_LISTE

pp = pprint.PrettyPrinter(indent=2)
import subprocess
import os
from subprocess import DEVNULL, STDOUT
#import pandas as pd
import numpy as np
from datetime import datetime as dt
import re
import time


class Phone:
    def __init__(self, id) -> None:
        self.id = id
        self.ip = f"{IP_BASE}.{id}"
        self.url = f"http://{self.ip}:{PHYPHOX_PORT}"
        print(f"[#{self.id}] Initializing... Available on network: {self.is_available}")

        self.port = 5400+id
        self.id_adb=ID_ADB_LISTE(id-100)


    @property
    def is_available(self):
        try:
            self.send_custom("config?", show_response=False, get_response=False)
            print(f"Phyphox connexion with phone {self.id} working")
            return True
        except:
            print(f"Phone {self.id} not phyphox connected")
            return False

    @property
    def is_running(self):
        req = self.send_custom("get?")["status"]
        is_running = req["measuring"]
        print(self.id)
        if is_running:
            if req["timedRun"]:
                print(f"Timed run -\t Remaining : {req['countDown']/1000:.2f} s")
            else:
                print("Running in non-stop mode...")
        else:
            print("Acquisition stopped")
        return is_running

    @property
    def start_time(self):
        return self.send_custom("time?=full")[0]["systemTime"]

    @property
    def is_adb(self):

        print(self.id_adb)
        c = subprocess.run(['adb','devices'],text=True,capture_output = True,shell=False)
        s = c.stdout.split('\n')
        lines = [line.split('\t') for line in s]
        if self.id_adb in lines:
            print(f"Phone {self.id} adb available")
        else:
            print()
        return None

    def launch_phyphox(self):
        subprocess.run(
            [
                "adb",
                "-s",
                f"{self.ip}:{ADB_PORT}",
                "shell",
                "monkey",
                "-p",
                "de.rwth_aachen.phyphox",
                "-c",
                "android.intent.category.LAUNCHER",
                "1",
            ],
            text=False,
            capture_output=False,
            stdout=DEVNULL,
            stderr=DEVNULL,
        )

    def get_start_time(self):
        self.start_time = self.send_custom("time?=full")[0]["systemTime"]

    def unlock(self):
        print("Unlocking phone")
        subprocess.run(
            ["bash", f"{BASH_SCRIPTS_PATH}/unlock.sh", f"{self.ip}:{ADB_PORT}"]
        )
        time.sleep(3)

    def activate_timedRun(self, t_before: int = 1, t_run: int = 100):
        """Runs a bash file to does action on the phoe to activate the Timed Run optio (Fixed Experiment duration)

        Args:
            - t_before (int, optional): Delay before the start of an experiment, in seconds. Defaults to 1.
            - t_run (int, optional): Experiment duration, in seconds. Defaults to 100.
        """
        print(
            f"Activating timed run with :\n\t- Delay before start:\t{t_before} s\n\t- Experiment duration:\t{t_run} s"
        )
        subprocess.run(
            [
                "bash",
                f"{BASH_SCRIPTS_PATH}/activate_timedRun.sh",
                f"{self.ip}:{ADB_PORT}",
                f"{t_before}",
                f"{t_run}",
            ]
        )
        time.sleep(10)

    def connect(self):
        print(f">> [#{self.id}] Initializing connexion...")
        ret = subprocess.run(
            ["adb", "connect", f"{self.ip}:{ADB_PORT}"],
            text=False,
            capture_output=False,
            stdout=DEVNULL,
            stderr=DEVNULL,
        )
        if ret != 0:
            print("Connexon successful ! ")

    def clear_buffers(self):
        with requests.get(f"{self.url}/control?cmd=clear"):
            print(f"\t\t[#{self.id}] >>  Clearing phone buffers")

    def start(self):
        with requests.get(f"{self.url}/control?cmd=start"):
            print(f"\t\t[#{self.id}] >>  Starting acquisition")

    def stop(self):
        with requests.get(f"{self.url}/control?cmd=stop"):
            print(f"\t\t[#{self.id}] >>  Stopping acquisition")

    def send_custom(
        self, request: str, show_response: bool = False, get_response: bool = True
    ):
        """Sends a customized request to a phone.

        Args:
            - request (str): the text of the request
            - show_response (bool, optional): Show in terminal the response obtained to the request. Defaults to False.
        """
        with requests.get(f"{self.url}/{request}") as response:
            if show_response:
                print(f'\t\t[#{self.id}] >>  Sending: "{request}"')
                print("Got response: ")
                pp.pprint(response.json())

            if get_response:
                return response.json()

    def dump(self, vars=["acc"], t=0):
        str_date = time.strftime("%d-%m-%Y>%H:%M:%S", time.gmtime(self.start_time))
        dims = ["X", "Y", "Z", "_time"]

        for var in vars:
            dimD_vars = [f"{var}{dim}" for dim in dims]
            cmd = "&&".join([f"{dimD_var}={t}|{var}_time" for dimD_var in dimD_vars])
            raw_container = self.send_custom(f"get?{cmd}")

            data_array = np.array(
                [raw_container["buffer"][dimD_var]["buffer"] for dimD_var in dimD_vars]
            ).T

            df = pd.DataFrame(data_array, columns=dimD_vars)

            df[f"{var}_time_abs"] = [
                dt.fromtimestamp(t + self.start_time) for t in df[f"{var}_time"]
            ]

            df.to_csv(os.path.join(DUMP, str_date, f"{str_date}_{var}_t={t:.2f}s.csv"))

    def start_datalog(self, var):

        self.stop()
        self.clear_buffers()
        self.start()

        time.sleep(1)

        str_date = time.strftime("%d-%m-%Y>%H:%M:%S", time.gmtime(self.start_time))
        os.makedirs(os.path.join(DUMP, str_date))
        while self.is_running:
            current_time = self.send_custom(f"get?{var}_time")["buffer"][f"{var}_time"][
                "buffer"
            ][0]
            self.dump(vars=[var], t=current_time - 1)

    def run_experiment(self, t_run):

        self.unlock()
        self.activate_timedRun(t_run=t_run)

        print(f"###########################\n[#{self.id}] >> Starting experiment\n")
        # exp = Experiment()

        self.clear_buffers()
        self.start()
        time.sleep(2)
        while self.is_running:
            pass

        self.datalog("acc")
