#!/usr/bin/env python3

import argparse
import os
import sys
import math
import csv


class Description:
    NAME = "discavar"
    VERSION = "0.1"
    SHORT = "Analyse multi-sample snp files with sample-relation " + \
            "information in vcf/vcf.gz format."


class tcolor:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def fmt_error(text):
    return(tcolor.FAIL +
           text +
           tcolor.ENDC)


def fmt_warning(text):
    return(tcolor.WARNING +
           text +
           tcolor.ENDC)


def fmt_success(text):
    return(tcolor.OKGREEN +
           text +
           tcolor.ENDC)


def fmt_info(text):
    return(text)


def warning(msg):
    print(fmt_warning("[WARNING] " + msg))


def error(msg):
    print(fmt_error("[ERROR] " + msg))
    sys.exit(1)


def success(msg):
    print(fmt_success("[OK] " + msg))


def info(msg):
    print(fmt_info("[INFO] " + msg))


class InputFile:
    def __init__(self,
                 filepath,
                 is_vcf=False,
                 is_tsv=False):
        self.filepath = os.path.abspath(filepath)
        self.exists = os.path.exists(self.filepath)
        self.is_vcf = self._is_vcf_file()
        self.is_tsv = self._is_tsv_file()

    def _is_vcf_file(self):
        vcf_endings = (".vcf", ".vcf.gz",
                       ".gvcf", ".gvcf.gz",
                       ".g.vcf", ".g.vcf.gz")
        if any(map(lambda e: self.filepath.endswith(e), vcf_endings)):
            return True
        return False

    def _is_tsv_file(self):
        tsv_endings = (".tsv")
        if any(map(lambda e: self.filepath.endswith(e), tsv_endings)):
            return True
        return False


class InputFileAction(argparse.Action):
    """
    Handles input files.
    Binds an InputFile object to the arparse namespace.
    """

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            if len(values) == 1:
                file = InputFile(values[0])
                if not file.exists:
                    msg = "File not found: {}".format(file.filepath)
                    raise(parser.error(msg))
                if not file.is_tsv and not file.is_vcf:
                    msg = "Invalid File: {}".format(file.filepath)
                    raise(parser.error(msg))
            else:
                msg = "More than one input file provided."
                raise(parser.error(msg))
        except RuntimeError:
            parser.error()

        # add the object to the namespace
        setattr(namespace, self.dest, file)


class OutdirAction(argparse.Action):
    """
    Handles outdirs mentionned as parameters
    """

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            if len(values) == 1:
                outdir = values[0]
                if os.path.exists(outdir):
                    subdir = "Discavar"
                    msg = "Outdir exists: {}".format(outdir) + \
                          "Using {} instead".format(subdir)
                    outdir = os.path.join(outdir, subdir)
            else:
                msg = "More than one output dir provided."
                raise(parser.error(msg))
        except RuntimeError:
            parser.error()

        # add to namespace
        setattr(namespace, self.dest, outdir)


class ChromRegionsAction(argparse.Action):
    """
    Handles interval files for filtering. Binds a list of ChromRegion
    objects to the argparse namespace
    """

    class ChromRegion:
        """
        Class Representing a Chromosomal Region.
        Initialized with:
        @param chrom    : Chromosome number
        @param start    : Start Position on that Chromosome. 
                          Default = 1
        @param end      : End Position on that chromosome
                          Default = math.infty
        """
        def __init__(self, chrom, start=1, end=math.inf):
            self.start = start
            self.end = end
            self.chrom = chrom

        def __str__(self):
            if self.start and self.end:
                region = "{}:{}-{}".format(self.chrom, self.start, self.end)
            else:
                region = self.chrom
            return region

    def __call__(self, parser, namespace, values, option_string=None):
        regions = []
        if values:

            if len(values) == 1 and self._is_gatk_file(values[0]):
                regions = self.read_gatk_file(values[0])

            elif len(values) == 1 and self._is_bed_file(values[0]):
                regions = self.read_bed_file(values[0])

            else:
                for region in values:
                    region = region.rstrip(",").rstrip(";").rstrip(" ")
                    regions.append(self.parse_gatk_str(region))

        # add the object to the namespace
        setattr(namespace, self.dest, regions)

    @classmethod
    def _is_gatk_file(self, filename):
        gatk_endings = (".list", ".intervals")
        if any(map(lambda e: filename.endswith(e), gatk_endings)):
            return True
        return False

    @classmethod
    def _is_bed_file(self, filename):
        bed_endings = (".bed")
        if any(map(lambda e: filename.endswith(e), bed_endings)):
            return True
        return False

    @classmethod
    def parse_gatk_str(self, gatk_region):
        """
        Takes a gatk style String describing a chromosomal region and returns
        a ChromRegion Object.
        @param gatk_region     : A String in the format <chr>:<start>-<end>
        """
        NOT_FOUND = -1
        region = None

        if gatk_region.find(":") != NOT_FOUND:
            chrom, pos = gatk_region.split(":")

            if pos.find("-") != NOT_FOUND:
                start, end = pos.split("-")
                start, end = int(start), int(end)
                region = self.ChromRegion(chrom, start, end)

        if not region:
            chrom = gatk_region
            region = self.ChromRegion(chrom=chrom)

        return region

    @classmethod
    def read_gatk_file(self, file):
        """
        Read all regions from gatk style file.
        Returns a list of ChromRegion objects.
        """
        regions = []
        with open(file, "r") as gatk_file:
            reader = csv.reader(gatk_file)
            for gatk_region in reader:
                regions.append(self.parse_gatk_str(gatk_region))

        return regions

    @classmethod
    def read_bed_file(self, file):
        """
        Read all chromosomal regions from a bed file.
        Returns a list of ChromRegion objects.
        """
        regions = []
        with open(file, "r") as bed_file:
            reader = csv.reader(bed_file, delimiter="\t")
            for line in reader:
                if line.startswith("#") or \
                   line.startswith("@"):
                    continue
                else:
                    if len(line) == 3:
                        chrom, start, end = line
                        region = self.ChromRegion(chrom, start, end)
                    else:
                        chrom = line[0]
                        region = self.ChromRegion(chrom)

                    regions.append(region)

        return region


def is_valid_job(job, parser):
    """
    Check if this job is sane.
    """

    if not job.subcommand:
        error_msg = "No analysis task selected. Has to be one of " + \
                    "vcfstats, cohort, versus"
        parser.error(error_msg)

    return True


def run_vcfstats():
    print("Running VCF statistics")


def parse_args(args=None):
    """
    Processes command line arguments.
    """

    parser = argparse.ArgumentParser(
            description=Description.SHORT,
            prog=Description.NAME)

    # version info
    parser.add_argument(
        "-V", "--version", action="version",
        version="%(prog)s "+Description.VERSION)

    subparsers = parser.add_subparsers(dest="subcommand")

    # parsing options
    # create db from snp files
    parser.add_argument(
        "--load-db", dest="create_db", action="store_true",
        help="Create db from snp files. SNP files will not be kept " +
        " in memory then. Recommended for large files/small ram.")

    # overview statistics
    vcfstat_parser = subparsers.add_parser(
        "vcfstats", help="VCF file statistics and filtering")
    # input options
    vcfstat_input = vcfstat_parser.add_argument_group("Input options")
    vcfstat_input.add_argument(
        "--vcf", type=str, action=InputFileAction, dest="input",
        default="", nargs=1, required=True,
        help="path to a single vcf/vcf.gz file")

    # familial analysis
    cohort_parser = subparsers.add_parser(
        "cohort", help="Familial/Cohort Analysis.")
    # input options
    cohort_input = cohort_parser.add_argument_group("Input options")
    cohort_input_me = cohort_input.add_mutually_exclusive_group(required=True)
    cohort_input_me.add_argument(
        "--vcf", type=str, action=InputFileAction, dest="input", nargs=1,
        help="path to a single vcf/vcf.gz file")
    cohort_input_me.add_argument(
        "--tsv", type=str, action=InputFileAction, dest="input", nargs=1,
        help="path to a tsv file containing paths to multiple" +
        "vcf/vcf.gz files and information on their relatedness")

    # healthy versus diseased analysis
    versus_parser = subparsers.add_parser(
        "versus", help="Healthy versus diseased analysis.")
    versus_input = versus_parser.add_argument_group("Input options")
    versus_input_me = versus_input.add_mutually_exclusive_group(required=True)
    versus_input_me.add_argument(
        "--vcf", type=str, nargs="?", dest="vcf_file",
        help="path to a single vcf/vcf.gz file")
    versus_input_me.add_argument(
        "--tsv", type=str, nargs=1, dest="tsv_file",
        help="path to a tsv file containing paths to multiple" +
        "vcf/vcf.gz files and information on their relatedness")

    # add filter options to every subparser
    for name, subp in subparsers.choices.items():

        # filter options
        filter_group = subp.add_argument_group("Variant filter options")
        filter_group.add_argument(
            "--quality-gt", type=float, dest="filter_gq",
            help="Filter variants below a quality score.")
        filter_group.add_argument(
            "--mapping-qual", type=float, dest="filter_mq",
            help="Filter variants below a given mapping quality")
        filter_group.add_argument(
            "--read-depth", type=float, dest="filter_dp",
            help="Filter variants with a read depth below the given threshold")
        filter_group.add_argument(
            "--regions", type=str, nargs="+",
            dest="filter_regions", default=[], action=ChromRegionsAction,
            help="Consider only variants in the defined region. Can be " +
            "a list of regions in the format 'chr2:340-100' or the path " +
            "to a file in GATK style or BED format.")
        filter_group.add_argument(
            "--ignore-regions", type=str, nargs="+",
            dest="filter_excl_regions", default=[], action=ChromRegionsAction,
            help="Ignore the variant in the given regions. Can be a " +
            "list of comma separated regions or the path to an interval " +
            "list file in GATK style or BED format")
        filter_group.add_argument(
            "--vtypes", dest="filter_vtypes", default=[],
            choices=["insertion", "deletion", "snp", "snv", "indel"],
            help="Consider only variants of a certain type. " +
            "For a explanation of variant types, please see: " +
            "http://vcftools.sourceforge.net/VCF-poster.pdf")
        filter_group.add_argument(
            "--ignore-types", dest="filter_excl_vtypes", default=[],
            choices=["insertion", "deletion", "snp", "snv", "indel"],
            help="Filter Variants of a certain type")
        filter_group.add_argument(
            "--perc-pass", dest="filter_perc_pass", type=float,
            default=0.9, help="Percentage of variants in a multi-sample " +
            "record that have to pass the filters.")

        # annotation filters
        ann_group = subp.add_argument_group("Annotation Filter Options")
        ann_group.add_argument(
            "--vep-consequences", type=str, nargs="+", dest="vep_consequence",
            default=[],
            help="Filter for variants based on VEP Consequence Annotation.")
        ann_group.add_argument(
            "--vep-impact", type=str, nargs="+", dest="vep_impact", default=[],
            help="Filter for variants based on VEP IMPACT Annotation")
        ann_group.add_argument(
            "--vep-significance", type=str, nargs="+", dest="vep_significance",
            default=[],
            help="Filter for variants based on VEP Significance annotation")

        # intersection options
        intersection_group = subp.add_argument_group("Intersection and " +
                                                     "subtractionOptions")
        intersection_group.add_argument(
            "--isec-alt-ratio", type=float, dest="alt_ratio", default=0.9,
            help="Alternate allel ratio for computing intersections")
        intersection_group.add_argument(
            "--isec-min-cr", type=float, dest="call_rate", default=0.9,
            help="Minimum callrate for a record to be considered " +
                 "an intersection")
        intersection_group.add_argument(
            "--sub-alt-ratio", type=float, dest="alt_ratio2", default=0.1,
            help="Maximum alternate allel ratio for the subtrahend group " +
            "when computing subtraction")
        intersection_group.add_argument(
            "--sub-max-cr", type=float, dest="call_rate2", default=0.5,
            help="Maximum callrate for the subtrahend group when computing " +
            "subtraction")

        # Output options
        output_group = subp.add_argument_group("Output options")
        output_group.add_argument(
            "-r", "--report-outdir", type=str, nargs=1,
            dest="report_dir", default="DiscavarReports", action=OutdirAction,
            help="Output directory for Reports.")
        output_group.add_argument(
            "-o", "--vcf-outdir", type=str, nargs=1,
            dest="vcf_dir", default="Discavar", action=OutdirAction,
            help="Output Directory for VCF files.")
        output_group.add_argument(
            "--vcf-copy", type=str, nargs=1,
            dest="vcf_copy", default="",
            help="Do not modify the input VCF and operate on a copy " +
                 "instead. The copy will be placed in the VCF outdir.")
        output_group.add_argument(
            "--interactive", action="store_true", dest="interactive_report",
            default=False,
            help="Export a JuPyter notebook with filtered vcf data. " +
            "NOT IMPLEMENTED")

    if args:
        job = parser.parse_args(args)

    else:
        job = parser.parse_args()

    if is_valid_job(job, parser):
        return job
