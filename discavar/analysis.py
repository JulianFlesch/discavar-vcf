from .vcfwrapper import VCFWrapper, get_samples_from_vcf
from .annotation_parser import VEPAnnotation
from .cli import InputFile, error, warning, info
from .filters import VariantFilter, AnnotationFilter

from cyvcf2 import VCF
from scipy import stats, mean
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
import allel
import re
import math
import os
import csv
import time


# static factory method for returning a vcf object
def load_vcf(vcf_file):
    pass
    if os.path.exists(vcf_file):
        vcf = VCF(vcf_file)
        return vcf
    else:
        return None


# static factory method for returning a cohort object representing
# family relations
def load_group(tsv_file):
    pass


def order_chroms(chroms):
    """
    static method for ordering a list of chromosome names correctly without
    knowledge of the chromosome name format (e.g. "chr1", "chr-1", or "CHROM1")

    @param chroms   : Iterable returning chromosome identifiers used in
                      a vcf file.

    @return ordered List of chromosomes
    """
    gonosomes = ["X", "Y"]
    mit_chrom = ["M", "MT"]
    ordered = []
    # find the chromosome-id prefix: i.e.the part that is shared among
    # all chromosome-ids
    chrom_prefix_re = re.compile(r"^[a-zA-Z_:\*\+\-]*")
    chrom_prefix = chrom_prefix_re.match(chroms[0]).group(0)
    # order in which the chromosomes will be printed:
    # first chr1 - chr23, then chrX and finally chrY
    ordered_chromnum = list(str(i) for i in range(1, 25)) \
        + gonosomes \
        + mit_chrom
    for i in ordered_chromnum:
        chrom = chrom_prefix + i
        if chrom in chroms:
            ordered.append(chrom)
        else:
            chrom = chrom_prefix + i.zfill(2)
            if chrom in chroms:
                ordered.append(chrom)
            else:
                # chromosome was not originally in the list -> skip it
                continue
    return ordered


def get_color_map(n, name="hsv"):
    """
    static method for generating a color map. This can be used for plots
    to create distinct colors.

    @param n    : Number of different colors to be returned

    @return a color map of distinct colors
    """
    return plt.cm.get_cmap(name, n)


class AnalysisControl:
    """
    Class for setting up and building cohorts for running analyses from
    command line input.

    @param job      : An argparse namespace obtained from cli.parse_args()
                      representing command line options controlling an
                      analysis task.
    @param analysis : One of ["versus", "cohort", "vcfstats"] representing
                      the Anaylsis taks to be run.
    """

    def __init__(self,
                 analysis,
                 job):

        # analysis task to be run
        self.analysis = analysis
        if self.analysis == "vcfstats":
            self.run = self.run_vcfstats
        elif analysis == "cohort":
            self.run = self.run_cohort_analysis
        elif analysis == "versus":
            self.run = self.run_versus_analysis
        else:
            error("Unexpected Subcommand: {}".format(analysis))

        # output options
        self.base_dir = os.getcwd()
        self.interactive_report = job.interactive_report
        self.report_dir = "DiscavarReports"
        if job.report_dir:
            self.report_dir = job.report_dir

        if not os.path.exists(self.report_dir):
            os.mkdir(self.report_dir)
        else:
            warning("Report outdir {} exists".format(self.report_dir))

        if not self.analysis == "vcfstats":
            self.vcf_dir = "Discavar"
            if job.vcf_dir:
                self.vcf_dir = job.vcf_dir

            if not os.path.exists(self.vcf_dir):
                os.mkdir(self.vcf_dir)
            else:
                warning("VCF outdir {} exists".format(self.vcf_dir))

            os.chdir(self.vcf_dir)

        # isec and subt options
        self.call_rate = job.call_rate
        self.call_rate2 = job.call_rate2
        self.alt_ratio = job.alt_ratio
        self.alt_ratio2 = job.alt_ratio2

        self.cohorts = []
        if job.input.is_tsv:
            self.cohorts = Cohort.from_tsv(job.input.filepath,
                                           build_cohort=True)

        elif job.input.is_vcf:

            self.cohorts = [Cohort.from_vcf("discavar-cohort",
                                            job.input.filepath,
                                            outfile=job.vcf_copy)]

        else:
            error("Invalid Input file ending.")

        if self.cohorts:

            # only apply filters when cohort analysis is selected
            if self.analysis == "cohort":
                info("Applying Filters")

                # set up filters
                var_filter = VariantFilter(
                    min_gq=job.filter_gq,
                    min_mq=job.filter_mq,
                    min_dp=job.filter_dp,
                    intervals=job.filter_regions,
                    excl_intervals=job.filter_excl_regions,
                    vtypes=job.filter_vtypes,
                    excl_vtypes=job.filter_excl_vtypes,
                    perc_pass=job.filter_perc_pass)

                # set up vep filter
                file = self.cohorts[0].vcf_file
                vcf = VCF(file)

                vep_parser = VEPAnnotation(list(vcf.header_iter()))
                vep_filter = AnnotationFilter(
                                vep_parser,
                                IMPACT__in=job.vep_impact,
                                CLIN_SIG__in=job.vep_significance)

                for c in self.cohorts:
                    c.add_filter(var_filter)
                    c.add_filter(vep_filter)
                    c.apply_filters()

        else:
            warning("No Cohorts or VCF specified.")

    def run_vcfstats(self):
        info("Gathering VCF Statistics")
        for c in self.cohorts:
            # Reporting
            info("Writing reports for {}".format(c.cohort_id))
            os.chdir(self.base_dir)
            os.chdir(self.report_dir)

            if not os.path.exists(c.cohort_id):
                os.mkdir(c.cohort_id)

            os.chdir(c.cohort_id)

            c.write(write_vcf=False,
                    write_reports=True)

    def run_cohort_analysis(self):
        info("Starting Cohort Analysis")
        for c in self.cohorts:
            os.chdir(self.base_dir)

            # Variant analysis
            info("Running cohort variant operations on {}".format(c.cohort_id))
            os.chdir(self.vcf_dir)
            affected = c.get_affected_sample_ids()
            unaffected = c.get_unaffected_sample_ids()

            time_start = time.time()
            if affected and unaffected:
                c.healthy_vs_diseased(self.alt_ratio,
                                      self.call_rate,
                                      self.alt_ratio2,
                                      self.call_rate2)

            else:
                c.intersect(self.alt_ratio, self.call_rate)

            time_stop = time.time()
            duration = (time_stop - time_start) / 60
            info("cohort analysis of {} was completed in {} minutes".format(
                 c.cohort_id, duration))

            c.write(write_vcf=True,
                    write_reports=False)

            # Reporting
            info("Writing reports for {}".format(c.cohort_id))
            os.chdir(self.base_dir)
            os.chdir(self.report_dir)

            if not os.path.exists(c.cohort_id):
                os.mkdir(c.cohort_id)

            os.chdir(c.cohort_id)

            c.write(write_vcf=False,
                    write_reports=True)

    def run_versus_analysis(self):
        #info("Starting Versus Analysis")
        error("Versus Analysis not implemented")


class Report:
    """
    Class for writing a report textfile with statistics gathered from a
    cohort object.

    @param cohort   : A Cohort object for which a report is written
    """
    def __init__(self,
                 cohort,
                 commandline="",
                 report_file=""):

        self.cohort = cohort
        self.report_file = report_file
        if not report_file:
            self.report_file = self.cohort.cohort_id + "_report.txt"

        self.report = ""

    def add_line(self, s="", line_prepend=""):
        """
        Adds a line to the report.
        @param s    : a string representing the line that should be added.
        @param line_prepend     : a string representing an arbitrary set of
                                  characters that is written before every line
        """
        s = line_prepend + s.replace("\n", "\n"+line_prepend)
        line = "{}\n".format(s)
        self.report += line

    def write(self):
        """
        Write the report to a plain text file.
        """
        with open(self.report_file, "w") as outfile:
            outfile.write(self.report)


class CohortReport(Report):
    def gather(self):
        self.add_line("Cohort Report for {}".format(self.cohort.cohort_id))
        self.add_line("{} Samples in this Cohort: {}".format(
                        len(self.cohort.get_sample_ids()),
                        self.cohort.get_sample_ids())
                      )
        self.add_line()

        # output gene List
        genes = self.cohort.gene_list()

        # Plotting
        data = self.cohort.to_dataframe()

        if data is not None and \
           isinstance(data, pd.DataFrame) and \
           not data.empty:
            _ = self.cohort.manhattan(df=data)
            _ = self.cohort.callrate_vs_qual(df=data)
            _ = self.cohort.qual_dist()
            dp_dist = self.cohort.dp_dist()
            mq_dist = self.cohort.mq_dist()
            vep_csq = self.cohort.vep_consequences()
            vtypes = self.cohort.vtypes_dist()
            # Todo: Order by Chromosomes
            var_stats = data.describe()
            #dp_stats = dp_dist.describe()
            #mq_stats = mq_dist.describe()

        else:
            error("No Statistics on the Variants could be collected")

        self.add_line("Variants with Cohort Relevance: {}".format(len(data)))
        self.add_line("Callrates:")
        self.add_line("\tMean:\t\t{}".format(var_stats["Callrate"]["mean"]))
        self.add_line("\tMax:\t\t{}".format(var_stats["Callrate"]["max"]))
        self.add_line("\tStd:\t\t{}".format(var_stats["Callrate"]["std"]))
        self.add_line("\t25%:\t\t{}".format(var_stats["Callrate"]["25%"]))
        self.add_line("\t50%:\t\t{}".format(var_stats["Callrate"]["50%"]))
        self.add_line("\t75%:\t\t{}".format(var_stats["Callrate"]["75%"]))
        self.add_line("")
        self.add_line("Alternate Allele Counts:")
        self.add_line("\tMean:\t\t{}".format(var_stats["AltAlleles"]["mean"]))
        self.add_line("\tMax:\t\t{}".format(var_stats["AltAlleles"]["max"]))
        self.add_line("\tStd:\t\t{}".format(var_stats["AltAlleles"]["std"]))
        self.add_line("\t25%:\t\t{}".format(var_stats["AltAlleles"]["25%"]))
        self.add_line("\t50%:\t\t{}".format(var_stats["AltAlleles"]["50%"]))
        self.add_line("\t75%:\t\t{}".format(var_stats["AltAlleles"]["75%"]))
        self.add_line()
        self.add_line("Variant Quality:")
        self.add_line("\tMean:\t\t{}".format(var_stats["QUAL"]["mean"]))
        self.add_line("\tMax:\t\t{}".format(var_stats["QUAL"]["max"]))
        self.add_line("\tStd:\t\t{}".format(var_stats["QUAL"]["std"]))
        self.add_line("\t25%:\t\t{}".format(var_stats["QUAL"]["25%"]))
        self.add_line("\t50%:\t\t{}".format(var_stats["QUAL"]["50%"]))
        self.add_line("\t75%:\t\t{}".format(var_stats["QUAL"]["75%"]))
        #self.add_line()
        #self.add_line("Mapping Quality:")
        #self.add_line(str(mq_stats))
        #self.add_line()
        #self.add_line("Read Depth:")
        #self.add_line(str(dp_stats))
        self.add_line()
        self.add_line("Variant Scores")
        self.add_line("\tMean:\t\t{}".format(var_stats["SCORE"]["mean"]))
        self.add_line("\tMin:\t\t{}".format(min(data["SCORE"])))
        self.add_line("\tMax:\t\t{}".format(var_stats["SCORE"]["max"]))
        self.add_line("\tStd:\t\t{}".format(var_stats["SCORE"]["std"]))
        self.add_line("\t25%:\t\t{}".format(var_stats["SCORE"]["25%"]))
        self.add_line("\t50%:\t\t{}".format(var_stats["SCORE"]["50%"]))
        self.add_line("\t75%:\t\t{}".format(var_stats["SCORE"]["75%"]))
        self.add_line()
        self.add_line("Variant Types")
        self.add_line(str(vtypes), line_prepend="\t")
        # TODO: Summary statistics
        self.add_line("VEP Consequences")
        self.add_line(str(vep_csq), line_prepend="\t")
        self.add_line()
        self.add_line("Genes Affected")
        unknown_gene_rowname = "unknown"
        num_affected = len(genes[genes["GENES"] != unknown_gene_rowname])
        self.add_line("\tTotal number affected:\t\t{}".format(num_affected))
        self.add_line("\tMost affected Genes:")
        self.add_line(str(genes.head(n=10)), line_prepend="\t")


class VersusReport(Report):
    def gather(self):
        pass


class Sample:
    """
    Class to represent a sample with a id, a vcf file, sex and disease
    status.
    """
    def __init__(self,
                 sample_id,
                 vcf_file,
                 sex,
                 affected):
        self.sample_id = sample_id
        self.vcf_file = vcf_file
        self.sex = sex
        self.affected = affected


class Cohort:
    """
    Class to represent a related cohort of Variants from different samples
    e.g. a family or healthy/diseased study.
        input_file      :   cli.InputFile Representation of the user input
    """
    def __init__(self,
                 cohort_id,
                 input_files=None,
                 filters=[],
                 regions=[],
                 info_rules=[],
                 use_database=False,
                 extra_threads=0,
                 build_cohort=False):

        self.cohort_id = cohort_id
        self.vcf_file = self.cohort_id + ".vcf.gz"
        # and instance of a VCFWrapper
        self.variants = None
        self.samples = []
        self.filters = filters
        self.regions = regions
        self.info_rules = info_rules
        self.build_cohort = False
        self.use_database = use_database
        self.extra_threads = extra_threads

    @classmethod
    def from_vcf(cls,
                 cohort_id,
                 file,
                 outfile="",
                 use_database=False):
        """
        Instantiate a Cohort from a single multisample vcf file.
        Assumes the vcf_file to contain a single cohort.
            vcf_file    :   String representing the path to a vcf_file
        """

        if not isinstance(file, InputFile):
            file = InputFile(file)

        if file.exists and file.is_vcf:
            info("Loading VCF {}".format(file.filepath))

            # new Cohort without VCFWrapper and Samples
            cohort = cls(cohort_id=cohort_id)

            # if an outfile is given
            if outfile != "":
                cohort.vcf_file = outfile

            else:
                cohort.vcf_file = file.filepath

            # build VCFWrapper
            variants = VCFWrapper.wrapper_factory(
                                    files=[file.filepath],
                                    cohort_filename=cohort.vcf_file,
                                    use_database=use_database)

            # extract sample ids from wrapper object
            sample_ids = variants.get_sample_names()

            # add VCFWrapper and Sample objects to the new cohort
            cohort.variants = variants
            for sid in sample_ids:
                cohort.add_sample(Sample(sample_id=sid,
                                         vcf_file=file.filepath,
                                         sex="NA",
                                         affected=False))

            return cohort

        else:
            warning("Invalid VCF Input: " +
                    "File doesn't exist or not in VCF Format")
            return None

    @classmethod
    def from_tsv(cls,
                 file,
                 build_cohort=True,
                 extra_threads=0,
                 use_database=False):
        """
        Instantiate a Cohort from a tsv file.
            tsv_file    :   String representing the path to a tsv file.
        """

        if not isinstance(file, InputFile):
            input_file = InputFile(file)

        cohorts = []

        if input_file.exists and input_file.is_tsv:
            info("Loading Cohorts from TSV file {}"
                 .format(input_file.filepath))

            with open(input_file.filepath, "r") as tsvfile:
                reader = csv.reader(tsvfile, delimiter="\t")

                # keep track of all  used sample ids for duplicate check
                used_sids = []

                # fetch all data from the tsv
                for line in reader:

                    # skip header lines
                    if line[0].startswith("#"):
                        continue

                    # unpack the line
                    sample_cohort = line[0]
                    sample_vcffile = os.path.abspath(line[3])
                    sample_status = line[1]
                    sample_sex = line[2]

                    if len(line) == 5:
                        # read sample ids
                        # the fifth column in the tsv is a comma separated list
                        # of ids.
                        sample_ids = [s.rstrip(" ").strip(" ")
                                      for s in line[4].split(",")]
                    else:
                        # get Sample id from vcf file
                        if os.path.exists(sample_vcffile):
                            sample_ids = get_samples_from_vcf(sample_vcffile)
                        else:
                            error("VCF file listed in {} not found: {}"
                                  .format(input_file.filepath,
                                          sample_vcffile))
                    while sample_ids:
                        sid = sample_ids.pop()

                        unique_sid = cls.get_unique_sid(
                                        sid, used_sids)

                        # reserve the sample_id that was constructed by adding
                        # it to the list of unique sample_ids for this cohort
                        used_sids.append(unique_sid)

                        # instatiate the sample
                        sample = Sample(vcf_file=sample_vcffile,
                                        sample_id=unique_sid,
                                        affected=sample_status,
                                        sex=sample_sex)

                        # if the cohort exists, add the sample object to it
                        for c in cohorts:
                            if c.cohort_id == sample_cohort:
                                c.add_sample(sample)
                                break

                        # else create a new cohort, add the sample to it,
                        # and add the new cohort to the list of cohorts
                        else:
                            cohort = cls(sample_cohort,
                                         build_cohort=build_cohort,
                                         extra_threads=extra_threads)
                            cohort.add_sample(sample)
                            cohorts.append(cohort)

        else:
            warning("Invalid TSV input: " +
                    "File doesn't exist or not in .tsv format")

        if build_cohort:
            if extra_threads:
                # build all cohorts at the same time
                pass
            else:
                for cohort in cohorts:
                    time_build_start = time.time()
                    cohort.build()
                    time_build_stop = time.time()
                    duration = (time_build_stop - time_build_start) / 60
                    info("Cohort building of {} completed in {} min".format(
                         cohort.cohort_id, duration))

        if cohorts:
            info("{} Samples were loaded into {} Cohorts"
                 .format(len(used_sids), len(cohorts)))
        else:
            warning("No Cohorts were build!")

        return cohorts

    @classmethod
    def get_unique_sid(cls, sid, used_sids):
        """
        Finds and returns a unique sample id
        """
        # find a unique sample ID for that cohort
        # If a sample id already exists inside a cohort,
        # a unique prefix of the form "<unique-int>:" is
        # prepended.
        # This mimicks the behavior of vcftools for dealing
        # with duplicate ids.
        unique_prefix = ""
        unique_sid = unique_prefix + sid

        # try all prefixes until one has not been used before
        while (unique_sid in used_sids):
            # construct a maybe unique prefix

            # A: No prefix used yet
            if unique_prefix == "":
                unique_prefix = "2:"

            # B: Or by incrementing the previous prefix
            else:
                unique_id = int(unique_prefix.strip(":")) + 1
                unique_prefix = str(unique_id) + ":"

            unique_sid = unique_prefix + sid

        if sid != unique_sid:
            warning("Sample ID '{}' was changed to the unique ID '{}'"
                    .format(sid, unique_sid))

        return unique_sid

    def get_sample_ids(self):
        """
        Return a list with the ids of all samples
        """
        return [s.sample_id for s in self.samples]

    def get_affected_sample_ids(self):
        """
        Returns a list with the ids of all affected samples
        """
        affected = []
        for s in self.samples:
            if s.affected:
                affected.append(s)

        return affected

    def get_unaffected_sample_ids(self):
        """
        Returns a list with the ids of all unaffected samples
        """
        unaffected = []
        for s in self.samples:
            if not s.affected:
                unaffected.append(s)

    def add_sample(self, sample):
        """
        Add sample object to representation of cohort
            sample      :   Sample representing a sample that should be added
        """
        self.samples.append(sample)

    def add_filter(self, fltr):
        """
        Add a filter of type Filter to the Cohort for building
        """
        self.filters.append(fltr)

    def add_filters(self, fltrs):
        """
        Add a set of filters to the Cohort for building
        """
        for fltr in fltrs:
            self.filters.append(fltr)

    def apply_filters(self, filters=[]):
        """
        Apply all filters.
        """
        if not filters:
            filters = self.filters

        time_start = time.time()

        self.variants.apply(filters)

        time_stop = time.time()
        duration = (time_stop - time_start) / 60
        info("Filtering of {} was completed in {} minutes".format(
             self.cohort_id, duration))

    def filter(self, filters=[]):
        """
        Query based on all filters.
        Returns an iterable of variants that passed the filter
        """
        return self.variants.query(filters=filters)

    def build(self):
        """
        Initialize the underlying wrapper object
        """
        filenames = [s.vcf_file for s in self.samples]
        self.variants = VCFWrapper.wrapper_factory(
                                    filenames,
                                    self.vcf_file,
                                    use_database=self.use_database)

        self.apply_filters(filters=self.filters)

    def intersect(self,
                  alt_ratio,
                  call_rate,
                  sample_groups=[]):
        """
        Intersect all samples in this cohort
        """
        if not sample_groups:
            # create one sample group containing all samples
            sample_groups = [
                [s.sample_id for s in self.samples]
            ]

        info("Computing intersection of {}".format(sample_groups))
        self.variants.intersect(sample_groups, alt_ratio, call_rate)

    def healthy_vs_diseased(self,
                            healthy_ar=1,
                            healthy_cr=1,
                            diseased_ar=0,
                            diseased_cr=0.5):
        """
        Perform a healthy versus diseased analysis of the samples in this
        cohort.
        """
        healthy = list(filter(lambda s: s.affected, self.samples))
        diseased = list(filter(lambda s: not s.affected, self.samples))

        self.variants.subtract(minuend=diseased,
                               subtrahend=healthy,
                               call_rate1=healthy_cr,
                               alt_ratio1=healthy_ar,
                               call_rate2=diseased_cr,
                               alt_ratio2=diseased_ar)

    def versus(self,
               cohort):
        """
        Versus Analysis of two cohort objects
        """
        pass

    def gene_list(self,
                  outfile=""):
        """
        Collects a list of genes and the number of variations that occur
        in each gene.
        Writes the List to a file
        """

        # currently, VEP Annotations are used to determine the gene name
        parser = self.variants.vep_parser
        if not parser:
            warning("Gene List can't be build without VEP Annotations")
            return None

        else:
            genes = {}

            for record in self.variants:

                # get the annotation string from the record
                ann_string = self.variants. \
                                get_var_info(record,
                                             field=parser.INFO_COL_NAME)

                # parse the string to obtain a dictionary representation
                # for every variant in the record as a list
                csq_annotations = parser.parse_annotation(ann_string)

                # it is assumed that the gene name is the same for all
                # variants, given that they all occur at the same 
                # chromosome and position
                csq_ann = csq_annotations[0]
                gene_name = csq_ann.get("Gene", "unknown")

                count = genes.get(gene_name, 0) + 1
                genes[gene_name] = count

            df = pd.DataFrame(dict(GENES=list(genes.keys()),
                                   VARIANTS=list(genes.values())))
            df.sort_values(by="VARIANTS", inplace=True, ascending=False)

            if not outfile:
                outfile = self.cohort_id + "_genelist.tsv"

            df.to_csv(outfile, sep="\t")

            return df

    def to_dataframe(self):
        """
        Build a Dataframe from the underlying vcf.
        This is more useful for plotting
        """
        fields = ["variants/REF",
                  "variants/ALT",
                  "variants/CHROM",
                  "variants/POS",
                  "variants/QUAL",
                  "calldata/GT"]

        data = allel.read_vcf(self.variants.cohort_filename,
                              fields=fields)

        if data is None or \
           not isinstance(data, dict):

            warning("Missing Variant data: " +
                    "No or incomplete data could be read from {}"
                    .format(self.variants.cohort_filename))

            return None

        try:
            df = pd.DataFrame(dict(CHROM=data["variants/CHROM"],
                                   POS=data["variants/POS"],
                                   QUAL=data["variants/QUAL"]))

            g = allel.GenotypeArray(data["calldata/GT"])

        except KeyError as e:
            print(e)
            warning("Missing Variant data: " +
                    "Not all required data could be read from {}"
                    .format(self.variants.cohort_filename))

        gn = g.to_n_alt()

        df["AltAlleles"] = [sum(vars) for vars in gn]
        df["Callrate"] = [(g.n_samples - list(vars).count(0)) / g.n_samples
                          for vars in gn]

        # compute variant score as the log product of their Quality and Callrate
        df["SCORE"] = [math.log2(qual * cr) for qual, cr in
                       zip(df["QUAL"], df["Callrate"])]

        # sort by variant position
        df.sort_values(["POS"])
        # sort again by chromosome name
        df["CHROM"] = df["CHROM"].astype("category")
        chromnames_ordered = order_chroms(df["CHROM"].unique())

        df["CHROM"] = df["CHROM"] \
            .cat.reorder_categories(chromnames_ordered)

        # add an index column
        df["IDX"] = range(len(df))

        return df

    def manhattan(self,
                  df=None,
                  outfile="",
                  only_data=False,
                  ax=None,
                  save_figure=True):
        """
        Manhattan plot or dataframe for the cohort.
        """

        if df is None or \
           type(df) != pd.DataFrame or \
           df.empty:
            df = self.to_dataframe()

        df_grouped = df.groupby("CHROM")

        if not only_data:
            # these colors will be repeated
            colors = get_color_map(23)

            if not ax:
                fig = plt.figure()
                ax = fig.add_subplot(111)

            else:
                fig = ax.get_figure()
                save_figure = False

            x_labels = []
            x_labels_pos = []
            for i, (name, group) in enumerate(df_grouped):
                if not "Y" in str(name).upper():
                    group.plot(kind='scatter',
                               x='IDX',
                               y='SCORE',
                               color=colors(i),
                               ax=ax)
                    x_labels.append(name)
                    x_labels_pos.append((group.IDX.iloc[-1] -
                                        (group.IDX.iloc[-1] -
                                         group.IDX.iloc[0]) / 2))

            ax.set_xticks(x_labels_pos)
            ax.set_xticklabels(x_labels)
            # roate the ticks, so each one is nicely readable
            for label in ax.get_xticklabels():
                label.set_rotation(50)
            ax.set_xlim([0, len(df)])
            ax.set_xlabel("Chromosome")
            ax.set_ylabel("Score")
            ax.set_title("Manhattan Distribution")

            # widen the figure
            fig.set_size_inches(15, 7)
            # fix cutoff labels
            fig.tight_layout()

            if not outfile:
                outfile = self.cohort_id + "_manhattan.png"

            if save_figure:
                fig.savefig(outfile)

        return df

    def vtypes_dist(self,
                    outfile="",
                    only_data=False,
                    ax=None,
                    save_figure=True):
        """
        Manhattan plot.
        """

        vtypes_hist = dict()

        # chromosomes listed in the vcf
        chroms = set()

        # vep parser for parse the VEP Annotation in the INFO Column
        parser = self.variants.vep_parser
        vep_annotated = True if parser else False

        for var in self.variants:

            var_chrom = self.variants.get_var_chrom(var)

            # save the chromosome the record references for plotting
            chroms.add(var_chrom)

            if vep_annotated:

                # read the raw VEP annotation string from the record INFO
                ann_string = self.variants. \
                                get_var_info(var,
                                             field=parser.INFO_COL_NAME)

                # parse the raw annotation string to obtain a dictionary
                # of <vep-annotation-field>:<vep-annotation> for each
                # variant in the record as a list.
                record_annotations = parser.parse_annotation(ann_string)

                for var_annotation in record_annotations:

                    var_type = var_annotation.get("VARIANT_CLASS", "unknown")

                    # count the variants number
                    chrom_hist = vtypes_hist.get(var_type, dict())
                    count = chrom_hist.get(var_chrom, 0) + 1
                    chrom_hist[var_chrom] = count
                    vtypes_hist[var_type] = chrom_hist

            # read variant type from the VCFWrapper
            else:
                # supports only one consensus variant type
                var_type = self.variants.get_var_type(var)

                # count the variants number
                chrom_hist = vtypes_hist.get(var_type, dict())
                count = chrom_hist.get(var_chrom, 0)
                chrom_hist[var_chrom] = count + 1
                vtypes_hist[var_type] = chrom_hist

        chroms = order_chroms(list(chroms))

        if not only_data:

            # plotting
            if not outfile:
                outfile = self.cohort_id + "_vartypes.png"
            if not ax:
                fig = plt.figure()
                ax = fig.add_subplot(111)
            else:
                fig = ax.get_figure()
                save_figure = False

            idx = range(len(chroms))
            bar_width = 0.8
            bar_below = []
            labels = []
            for vtype in vtypes_hist.keys():
                labels.append(vtype)
                freqs = [vtypes_hist[vtype].get(c, 0) for c in chroms]
                if bar_below:
                    ax.bar(idx,
                           freqs,
                           width=bar_width,
                           bottom=bar_below)
                    bar_below = [sum(x) for x in zip(freqs, bar_below)]
                else:
                    ax.bar(idx, freqs, width=bar_width)
                    bar_below = freqs

            ax.set_ylabel("Frequency")
            ax.set_xlabel("Chromosome")
            ax.set_xticklabels(chroms, rotation=50)

            ax.legend(labels)

            ax.set_title("Variant Types by Chromosome")

            # fix cutoff labels
            fig.tight_layout()

            if save_figure:
                fig.savefig(outfile)

        return pd.DataFrame(vtypes_hist).fillna(0)

    def callrate_vs_qual(self,
                         df=None,
                         outfile="",
                         ax=None,
                         save_figure=True):

        if df is None or \
           type(df) != pd.DataFrame or \
           df.empty:
            df = self.to_dataframe()

        df_grouped = df.groupby("CHROM")
        colors = get_color_map(len(df_grouped))
        legend = []

        if not outfile:
            outfile = self.cohort_id + "_callrate_v_vqual.png"

        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)

        else:
            fig = ax.get_figure()
            save_figure = False

        for i, (name, group) in enumerate(df_grouped):
            group.plot(kind="scatter",
                       x="Callrate",
                       y="QUAL",
                       color=colors(i),
                       ax=ax)
            legend.append(name)

        ax.legend(legend, ncol=3)
        ax.set_xlabel("Callrate")
        ax.set_ylabel("Record mean QUAL")
        ax.set_ymargin(0.3)
        ax.set_title("{} variant relevance and mean quality"
                     .format(self.cohort_id))

        # fix cutoff labels
        fig.tight_layout()

        if save_figure:
            fig.savefig(outfile)

        return df

    def dp_dist(self,
                outfile="",
                ax=None,
                save_figure=True,
                bins=50):

        DP_field = "DP"
        dps = []
        for var in self.variants:
            dps.append(self.variants.get_var_info(var, field=DP_field))

        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title("Read Depth Distribution")

        else:
            fig = ax.get_figure()
            save_figure = False

        ax.set_autoscaley_on(True)
        ax.hist(dps, bins=bins)

        # logarithmic y scale
        ax.set_yscale("log")

        ax.set_xlabel("Read Depth")
        ax.set_ylabel("Frequency")
        ax.set_title("Read Depth Distribution")

        # fix cutoff labels
        fig.tight_layout()

        if not outfile:
            outfile = self.cohort_id + "_dp_dist.png"

        if save_figure:
            fig.savefig(outfile)

        return dps

    def mq_dist(self,
                outfile="",
                ax=None,
                save_figure=True,
                bins=50):

        MQ_field = "MQ"
        mqs = []
        for var in self.variants:
            mqs.append(self.variants.get_var_info(var, field=MQ_field))

        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title("Mapping Quality Distribution")

        else:
            fig = ax.get_figure()
            save_figure = False

        ax.hist(mqs, bins=bins)

        # logarithmic y scale
        ax.set_yscale("log")

        ax.set_xlabel("Mapping Quality")
        ax.set_ylabel("Frequency")
        ax.set_title("Mapping Quality Distribution")

        # fix cutoff labels
        fig.tight_layout()

        if not outfile:
            outfile = self.cohort_id + "_mq_dist.png"

        if save_figure:
            fig.savefig(outfile)

        return mqs

    def qual_dist(self,
                  outfile="",
                  ax=None,
                  save_figure=True,
                  bins=50):
        quals = []
        for var in self.variants:
            quals.append(mean(self.variants.get_var_gq(var)))

        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111)
        else:
            fig = ax.get_figure()
            save_figure = False

        ax.hist(quals, bins=bins)

        # logarithmic y scale
        ax.set_yscale("log")

        ax.set_xlabel("Quality")
        ax.set_ylabel("Frequency")
        ax.set_title("Genotype Quality Distribution")

        # fix cutoff labels
        fig.tight_layout()

        if not outfile:
            outfile = self.cohort_id + "_qual_dist.png"

        if save_figure:
            fig.savefig(outfile)

        return quals

    def vep_consequences(self,
                         show=8,
                         outfile="",
                         ax=None,
                         save_figure=True):
        """
        Plot the VEP Consequence annotations as a pie chart.
        Requires VEP Annotations.
        @param show     : integer describing how many of the most frequent
                          consequences should be shown.

        """

        if not self.variants.vep_parser:
            warning("VEP Consequence Plot requires VEP Annotation")
            return None

        else:
            consequences = {}

            # reference the vep_parser in self.variants for better
            # readability
            parser = self.variants.vep_parser
            for record in self.variants:

                # get the raw annotation string for the field "CSQ" in INFO for
                # a record
                ann_string = self.variants.get_var_info(record, field="CSQ")

                # obtain a dictionary representation of
                # <vep-fieldname>: <vep-annotation-value> pairs for every
                # variant in the record as a list
                record_annotations = parser.parse_annotation(ann_string)

                # here, the consequences for each variant in a record
                # are considered separately, since they might differ.
                for variant_ann in record_annotations:

                    cons = variant_ann.get("Consequence", "unknown")

                    # if more than one consequence is reported for a variant,
                    # they are returned as a list
                    if isinstance(cons, list):

                        for c in cons:
                            # retrieve the number of time this consequence
                            # has been called and increment by one
                            count = consequences.get(c, 0) + 1
                            # write again to the dict with results
                            consequences[c] = count

                    else:

                        # proceed like above except for single consequence
                        count = consequences.get(cons, 0) + 1
                        consequences[cons] = count

            df = pd.DataFrame.from_dict(consequences,
                                        orient="index",
                                        columns=["count"])

            # sort dataframe by count
            df = df.sort_values("count", ascending=False)

            if not ax:
                fig = plt.figure()
                ax = fig.add_subplot(211)
            else:
                fig = ax.get_figure()
                save_figure = False

            # elect the biggest categories
            most_freq_csq = df.head(show)

            num_other = df["count"].sum() - \
                most_freq_csq["count"].sum()

            most_freq_csq = most_freq_csq.append(
                pd.DataFrame.from_dict({"other": {"count": num_other}},
                                       orient="index"))

            labels = list(most_freq_csq.index)

            # plotting
            ax = most_freq_csq.plot.barh()

            ax.set_title("VEP Consequences")
            ax.set_xlabel("Frequency")
            ax.get_legend().remove()

            # change x ticks
            x_ticks = ticker.ScalarFormatter(useOffset=True)
            x_ticks.set_powerlimits((-3, 4))
            ax.xaxis.set_major_formatter(x_ticks)

            #ax.set_xticklabels(labels, rotation=50)

            fig = ax.get_figure()
            fig.tight_layout()
            if save_figure:
                if not outfile:
                    outfile = self.cohort_id + "_vep_csq.png"

                fig.savefig(outfile)

            return df

    def pca_cluster(self,
                    outfile=""):
        """
        Perform PCA on
        """
        fields = ["calldata/GT"]
        data = allel.read_vcf(self.variants.cohort_filename, fields=fields)
        g = allel.GentotypeChunkedArray(data["calldata/GT"])
        gn = g.to_n_alt()
        gn_clean = gn[~np.isnan(gn)]
        coords, model = allel.pca(gn_clean,
                                  n_components=4,
                                  scaler="patterson")

        fig = plt.figure()
        # pc1 vs pc2
        ax = fig.add_subplot(121)
        ax.scatter(coords[:, 0], coords[:, 1])
        # pc3 vs pc4
        ax = fig.add_subplot(122)
        ax.scatter(coords[:, 2], coords[:, 3])

    def write(self, write_vcf=True, write_reports=True):
        """
        Writes Result files:
            Cohort VCF with all relevant Variants and a cohort annotation
            Optional Report Files:
                Report TSV file with variants/genes
        """
        if write_vcf:
            self.variants.to_vcf()

        if write_reports:
            time_start = time.time()

            report = CohortReport(self)
            report.gather()
            report.write()

            time_stop = time.time()
            duration = (time_stop - time_start) / 60
            info("Reporting for {} was completed in {} minutes".format(
                 self.cohort_id, duration))


class SingleVCFStatistics:
    """
    Class to perform single VCF single Sample coverage statistics for a
    quick overview of the coverage and spread of variants.

        vcf_object:     cyvcf2.VCF Representation of a vcf

    """
    VARIANT_TYPES = {"snp", "indel", "unknown"}

    def __init__(self, vcf_object):
        self.vcf = vcf_object
        self.samples = self.vcf.samples
        self.chromosomes = None
        self._extract()

    def _extract(self):
        chromosomes = {}
        for var in self.vcf:
            if var.CHROM in chromosomes.keys():
                chromosomes[var.CHROM][var.var_type] += 1
            else:
                chromosomes[var.CHROM] = dict(
                    [(var, 0) for var in self.VARIANT_TYPES])
                chromosomes[var.CHROM][var.var_type] = 1

        self.chromosomes = pd.DataFrame(chromosomes)
        self.chromosomes.fillna(0)

    def print_stats(self):
        print("VCF file statistics")
        print("Variation Table:")
        print("Chrom\tindels\tsnps\tunknown\t")
        for chrom in self.chromosomes:
            print("{}:\t{}\t{}\t{} ".format(
                chrom,
                self.chromosomes[chrom]["indel"],
                self.chromosomes[chrom]["snp"],
                self.chromosomes[chrom]["unknown"]))

        row_sums = self.chromosomes.sum(axis=1)

        print("SUM:\t{}\t{}\t{}".format(
            row_sums["indel"],
            row_sums["snp"],
            row_sums["unknown"]))

        print("Total Variations: {}".format(
            sum([row_sums["indel"],
                row_sums["snp"],
                row_sums["unknown"]])))
