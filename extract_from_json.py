"""Program for the json files from https://neuromancer.sk/std/other/ using main.py and tkinter (both not in standard lib).
Upon execution of this program, a file explorer should open where one can select a json file previously downloaded from the site or structured in the same way.
Some exceptions exist to the standard form, that can't be opened currently, e.g. has no generator specified.
"""

import json
import sys
import tkinter as tk
from tkinter import filedialog

import curves
import main as m


def main():
    try:
        file_object = get_file_object()

    except Exception:
        # Dealing with CWE-248, but note the problems of CWE-396, keeping it dev for the following file processing.
        print("File could not be opened.")
        sys.exit()

    data = json.load(file_object)

    curves_types = {
        'Weierstrass': curves.ShortWeierstrass,
        'Montgomery': curves.Montgomery,
        'Edwards': curves.Edwards,
        'TwistedEdwards': curves.TwistedEdwards
    }

    curve_type = data["form"]

    coe1 = curves_types[curve_type].coefficient1_name
    coe2 = curves_types[curve_type].coefficient2_name

    if curve_type == 'Montgomery':
        # Has parameters a,b in json, so not compliant with convention (and therefore curves.py).
        # Solution: Overwrite
        coe1 = 'a'
        coe2 = 'b'

    curve_arguments_dic = {
        "p": int(data["field"]["p"], 16),
        curves_types[curve_type].coefficient1_name: int(data["params"][coe1]["raw"], 16),
        curves_types[curve_type].coefficient2_name: int(data["params"][coe2]["raw"], 16),
        curves_types[curve_type].coordinate1_name: int(data["generator"]["x"]["raw"], 16),
        curves_types[curve_type].coordinate2_name: int(data["generator"]["y"]["raw"], 16)
    }

    user_curve = m.init_curve(curve_arguments_dic, curves_types, curve_type)

    m.ask_transtype_and_transform(user_curve)


def get_file_object():
    root = tk.Tk()
    root.withdraw()

    file_path_tuple = filedialog.askopenfilenames()

    file_path = file_path_tuple[0]

    print(f"You are opening {file_path}.")

    file_object = open(file_path)

    return file_object


if __name__ == "__main__":
    main()

