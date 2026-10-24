__prog__ = "initialise_funcs.py"
__version__ = "1.0.0"
# Author:      Michal Martin & Partha Talukdar
# Created:     09/06/2026
# Version history
# ---------------

import json
import os
from glob import glob
from time import sleep
import sys


sleepTime = 5


def initiation(form):
    """
    This function is called to initiate the programme to process non-GUI settings.
    """
    if not _read_setup_file(form):
        print("Setup failed")
        sleep(sleepTime)
        sys.exit(0)

    form.config_file = os.path.join(form.config_dir, "plink_utils_config.txt")


def _write_default_config_file(config_file):
    """
    Create a default config file if one does not already exist.
    """
    config_dir, _ = os.path.split(config_file)
    if config_dir and not os.path.isdir(config_dir):
        os.mkdir(config_dir)

    default_config = {
        "fnames": {
            "inp_fname": "",
            "out_dir": "",
            "root_name": "",
            "start_from_rec": "0",
            "num_recs": "99999999",
        }
    }

    with open(config_file, "w") as fconfig:
        json.dump(default_config, fconfig, indent=2, sort_keys=True)

    return default_config


def _read_setup_file(form):
    """
    Locate the PLINK executable, logo, configuration directory,
    and PED/MAP input files relative to this Python script.
    """

    print("USING SCRIPT-RELATIVE SETUP")

    script_dir = os.path.dirname(
        os.path.abspath(__file__)
    )

    plink_dir = os.path.join(
        script_dir,
        "plink"
    )

    form.plink_exe = os.path.join(
        plink_dir,
        "exe_dir",
        "plink2.exe"
    )

    form.fname_png = os.path.join(
        plink_dir,
        "Images",
        "plink_workflow_logo.png"
    )

    form.config_dir = os.path.join(
        plink_dir,
        "setup_config_files"
    )

    map_ped_dir = os.path.join(
        plink_dir,
        "map_ped"
    )

    required_paths = [
        (
            form.plink_exe,
            "PLINK executable"
        ),
        (
            form.fname_png,
            "GUI logo"
        ),
        (
            form.config_dir,
            "configuration directory"
        ),
        (
            map_ped_dir,
            "PED/MAP directory"
        ),
    ]

    for required_path, description in required_paths:
        if not os.path.exists(required_path):
            print(
                f"{description} was not found: "
                f"{required_path}"
            )
            return False

    ped_files = glob(
        os.path.join(map_ped_dir, "*.ped")
    )

    map_files = glob(
        os.path.join(map_ped_dir, "*.map")
    )

    if not ped_files or not map_files:
        print(
            "No matching PED and MAP files were found in:",
            map_ped_dir
        )
        return False

    ped_roots = {
        os.path.splitext(path)[0]
        for path in ped_files
    }

    map_roots = {
        os.path.splitext(path)[0]
        for path in map_files
    }

    matching_roots = sorted(
        ped_roots.intersection(map_roots)
    )

    if not matching_roots:
        print(
            "PED and MAP files with the same root "
            "filename were not found in:",
            map_ped_dir
        )
        return False

    if len(matching_roots) > 1:
        print(
            "Multiple PED/MAP datasets were found. "
            "Using:",
            matching_roots[0]
        )

    form.map_ped_file = matching_roots[0]

    return True


def write_config_file(form):
    """
    Write current selections to the config file.
    """
    config_file = form.config_file
    config = {
        "fnames": {
            "plink_chr": form.w_combo01.currentText(),
            "out_dir": form.w_lbl_outdir.text(),
            "bp_lower": form.w_bp_lower.text(),
            "bp_upper": form.w_bp_upper.text(),
            "root_name": form.w_root_name.text(),
        }
    }

    with open(config_file, "w") as fconfig:
        json.dump(config, fconfig, indent=2, sort_keys=True)


def read_config_file(form):
    """
    Read the last saved GUI selections from the config file, creating defaults
    if the file does not already exist.
    """
    config_file = form.config_file
    if not os.path.isfile(config_file):
        config = _write_default_config_file(config_file)
    else:
        with open(config_file, "r") as fconfig:
            config = json.load(fconfig)

    fnames = config.get("fnames", {})

    plink_chr = fnames.get("plink_chr")
    if plink_chr:
        idx = form.w_combo01.findText(str(plink_chr))
        if idx >= 0:
            form.w_combo01.setCurrentIndex(idx)

    form.w_lbl_outdir.setText(fnames.get("out_dir", ""))
    form.w_bp_lower.setText(fnames.get("bp_lower", ""))
    form.w_bp_upper.setText(fnames.get("bp_upper", ""))
    form.w_root_name.setText(fnames.get("root_name", ""))


def _write_default_setup_file(setup_file):
    """
    Stanza if setup_file needs to be created.
    """
    root_dir = os.getcwd()

    default_setup = {
        "setup": {
            "root_dir": root_dir,
            "config_dir": root_dir,
            "fname_png": root_dir,
            "plink_exe": root_dir,
            "map_ped": root_dir,
        }
    }

    with open(setup_file, "w") as fsetup:
        json.dump(default_setup, fsetup, indent=2, sort_keys=True)

    return default_setup

