__prog__ = "PlinkUtilsGUI.py"
__version__ = "1.1.0"
__author__ = "Michal Martin & Partha Talukdar"

import csv
import os
import sys
import threading
import time

import psutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if SCRIPT_DIR not in sys.path:
    sys.path.insert(0, SCRIPT_DIR)

from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from initialise_funcs import initiation, read_config_file, write_config_file
from plink_funcs import (
    extract_markers,
    formal_haplotype_block_detection,
    ld_decay,
    ld_decay_graph,
    potential_haplotype,
    snp_trait_analysis,
    two_snp_pair_anova,
    validate_inputs,
)

PERFORMANCE_SAMPLING_INTERVAL_SECONDS = 0.1


def _combined_resident_memory_mib(process):
    """Return simultaneous RSS for this Python process and its descendants."""
    processes = [process]
    try:
        processes.extend(process.children(recursive=True))
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        pass

    total_bytes = 0
    observed_pids = set()
    for item in processes:
        if item.pid in observed_pids:
            continue
        observed_pids.add(item.pid)
        try:
            total_bytes += item.memory_info().rss
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
    return total_bytes / (1024 ** 2)


def _measure_operation(operation_name, operation, form):
    """Measure wall time and sampled peak RSS while a GUI operation runs."""
    process = psutil.Process(os.getpid())
    stop_event = threading.Event()
    baseline_mib = _combined_resident_memory_mib(process)
    peak_mib = baseline_mib

    def sample_memory():
        nonlocal peak_mib
        while not stop_event.is_set():
            peak_mib = max(peak_mib, _combined_resident_memory_mib(process))
            stop_event.wait(PERFORMANCE_SAMPLING_INTERVAL_SECONDS)

    sampler = threading.Thread(target=sample_memory, daemon=True)
    start_time = time.perf_counter()
    status = "returned"
    sampler.start()

    try:
        return operation(form)
    except Exception:
        status = "raised exception"
        raise
    finally:
        elapsed_seconds = time.perf_counter() - start_time
        stop_event.set()
        sampler.join(timeout=1.0)
        peak_mib = max(peak_mib, _combined_resident_memory_mib(process))

        output_directory = form.w_lbl_outdir.text()
        measurement_file = os.path.join(
            output_directory, "performance_measurements.csv"
        )
        fieldnames = [
            "operation",
            "elapsed_seconds",
            "sampled_peak_rss_mib",
            "baseline_rss_mib",
            "sampling_interval_seconds",
            "status",
        ]
        row = {
            "operation": operation_name,
            "elapsed_seconds": "{:.6f}".format(elapsed_seconds),
            "sampled_peak_rss_mib": "{:.3f}".format(peak_mib),
            "baseline_rss_mib": "{:.3f}".format(baseline_mib),
            "sampling_interval_seconds": PERFORMANCE_SAMPLING_INTERVAL_SECONDS,
            "status": status,
        }

        file_exists = os.path.isfile(measurement_file)
        with open(measurement_file, "a", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            if not file_exists:
                writer.writeheader()
            writer.writerow(row)

        print("\nPERFORMANCE MEASUREMENT")
        print("Operation: {}".format(operation_name))
        print("Elapsed time: {:.3f} seconds".format(elapsed_seconds))
        print("Sampled peak resident memory: {:.1f} MiB".format(peak_mib))
        print("Baseline resident memory: {:.1f} MiB".format(baseline_mib))
        print(
            "Requested sampling interval: {:.1f} seconds".format(
                PERFORMANCE_SAMPLING_INTERVAL_SECONDS
            )
        )
        print("Status: {}; success requires output verification".format(status))
        print("Measurements saved to: {}".format(measurement_file))

class Form(QWidget):
    def __init__(self, parent=None):
        super(Form, self).__init__(parent)

        self.version = "Plink utilities"
        initiation(self)

        hbox = QHBoxLayout()
        hbox.setSpacing(10)

        lh_vbox = QVBoxLayout()

        lbl20 = QLabel()
        pixmap = QPixmap(self.fname_png)
        pixmap = pixmap.scaledToWidth(340, Qt.SmoothTransformation)
        lbl20.setPixmap(pixmap)
        lh_vbox.addWidget(lbl20)

        hbox.addLayout(lh_vbox)

        rh_vbox = QVBoxLayout()

        grid = QGridLayout()
        grid.setSpacing(10)

        lbl01 = QLabel("Chromosome")
        help_text = "Select a particular chromosome"
        lbl01.setToolTip(help_text)
        grid.addWidget(lbl01, 1, 0)

        combo01 = QComboBox()
        combo01.setMaxVisibleItems(12)
        for chr_num in range(1, 13):
            combo01.addItem(str(chr_num))
        combo01.setMinimumContentsLength(3)
        combo01.setFixedWidth(80)
        combo01.view().setMinimumWidth(80)
        combo01.setToolTip(help_text)
        grid.addWidget(combo01, 1, 1)
        self.w_combo01 = combo01

        w_lbl05a = QLabel("Base pair lower:")
        grid.addWidget(w_lbl05a, 3, 0)

        w_bp_lower = QLineEdit()
        w_bp_lower.setFixedWidth(90)
        grid.addWidget(w_bp_lower, 3, 1)
        self.w_bp_lower = w_bp_lower

        w_lbl05b = QLabel("Base pair upper")
        help_text = "Maximum upper SNP value is 43,270,923"
        w_lbl05b.setToolTip(help_text)
        grid.addWidget(w_lbl05b, 3, 2)

        w_bp_upper = QLineEdit()
        w_bp_upper.setFixedWidth(90)
        w_bp_upper.setToolTip(help_text)
        grid.addWidget(w_bp_upper, 3, 3)
        self.w_bp_upper = w_bp_upper

        w_spacer04 = QLabel("")
        grid.addWidget(w_spacer04, 4, 1)

        w_overwrite_res = QCheckBox("Overwrite previously written results")
        w_overwrite_res.setChecked(True)
        grid.addWidget(w_overwrite_res, 8, 0, 1, 3)
        self.w_overwrite_res = w_overwrite_res

        w_out_dir = QPushButton("Output Directory")
        help_text = "Directory into which output files will be written"
        w_out_dir.setToolTip(help_text)
        grid.addWidget(w_out_dir, 10, 0)
        w_out_dir.clicked.connect(self.fetchOutDir)

        w_lbl_outdir = QLabel("")
        grid.addWidget(w_lbl_outdir, 10, 1, 1, 6)
        self.w_lbl_outdir = w_lbl_outdir

        w_lbl12 = QLabel("Root filename")
        grid.addWidget(w_lbl12, 12, 0)

        w_root_name = QLineEdit()
        w_root_name.setFixedWidth(200)
        grid.addWidget(w_root_name, 12, 1, 1, 2)
        self.w_root_name = w_root_name

        w_spacer18 = QLabel("")
        grid.addWidget(w_spacer18, 17, 1)

        w_spacer19 = QLabel("")
        grid.addWidget(w_spacer19, 18, 1)

        w_extract_markers = QPushButton("PLINK Extract")
        w_extract_markers.setToolTip("Extract markers")
        grid.addWidget(w_extract_markers, 19, 0)
        w_extract_markers.clicked.connect(self.markerExtraction)

        w_ld_decay = QPushButton("LD Decay")
        w_ld_decay.setToolTip("Run LD decay analysis")
        grid.addWidget(w_ld_decay, 19, 1)
        w_ld_decay.clicked.connect(self.runLdDecay)

        w_ld_graph = QPushButton("LD Decay Graph")
        w_ld_graph.setToolTip("Create LD decay PNG graph from the .vcor file")
        grid.addWidget(w_ld_graph, 19, 2)
        w_ld_graph.clicked.connect(self.runLdDecayGraph)

        w_snp_trait = QPushButton("Exploratory SNP-Trait Analysis")
        w_snp_trait.setToolTip("Run homozygote-only regional SNP-trait ANOVA with completeness and group-size filters")
        grid.addWidget(w_snp_trait, 19, 3)
        w_snp_trait.clicked.connect(self.runSnpTraitAnalysis)

        w_haplotype = QPushButton("Candidate LD Regions")
        w_haplotype.setToolTip("Summarise spatially adjacent strong-LD pairs as candidate LD regions")
        grid.addWidget(w_haplotype, 20, 0)
        w_haplotype.clicked.connect(self.runPotentialHaplotype)

        w_fhbd = QPushButton("Graph-Based LD Clusters")
        w_fhbd.setToolTip("Group strong-LD SNP relationships into graph-based connected components")
        grid.addWidget(w_fhbd, 20, 1)
        w_fhbd.clicked.connect(self.runFHBD)

        w_two_snp = QPushButton("Pairwise Genotype Combinations")
        w_two_snp.setToolTip("Compare grain Fe among homozygous two-locus genotype-combination groups")
        grid.addWidget(w_two_snp, 20, 2)
        w_two_snp.clicked.connect(self.runTwoSnpAnova)

        w_exit = QPushButton("Exit", self)
        w_exit.clicked.connect(self.exitClicked)
        grid.addWidget(w_exit, 20, 3)

        rh_vbox.addLayout(grid)
        hbox.addLayout(rh_vbox)
        self.setLayout(hbox)

        self.setGeometry(300, 300, 1040, 340)
        self.setWindowTitle(self.version + " - Integrated Regional Analysis")

        read_config_file(self)

    def fetchOutDir(self):
        fname_cur = self.w_lbl_outdir.text()
        fname = QFileDialog.getExistingDirectory(self, "Select directory", fname_cur)
        if fname != "" and fname != fname_cur:
            fname = os.path.normpath(fname)
            self.w_lbl_outdir.setText(fname)

    def markerExtraction(self):
        if validate_inputs(self):
            extract_markers(self)

    def runLdDecay(self):
        if validate_inputs(self):
            ld_decay(self)

    def runLdDecayGraph(self):
        if validate_inputs(self):
            ld_decay_graph(self)

    def runPotentialHaplotype(self):
        if validate_inputs(self):
            potential_haplotype(self)

    def runFHBD(self):
        if validate_inputs(self):
            formal_haplotype_block_detection(self)

    def runTwoSnpAnova(self):
        if validate_inputs(self):
            _measure_operation(
                "pairwise_genotype_combination", two_snp_pair_anova, self
            )

    def runSnpTraitAnalysis(self):
        if validate_inputs(self):
            _measure_operation(
                "exploratory_regional_snp_trait", snp_trait_analysis, self
            )

    def exitClicked(self):
        write_config_file(self)
        self.close()


def main():
    app = QApplication(sys.argv)
    form = Form()
    form.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()




