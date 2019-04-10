from .cli import error
from .annotation_parser import VEPAnnotation

from cyvcf2 import VCF, Writer
# wormtable
# imort wormtable
import os
import subprocess
import shutil


# static method for dependency checking.
def check_dependency(dependency, dependency_name=None):
    """
    Throws error, if the dependency does not exist or is not readable.
    @param dependency   : string representing the name of a system dependency
                            that needs to be installed
    """
    if not dependency_name:
        dependency_name = dependency

    path_to_binary = shutil.which(dependency)

    if not (path_to_binary and os.access(path_to_binary, os.R_OK)):
        error("Required dependency not satisfied: {}".format(dependency_name))


# static method for adding an extra extension to vcf filenames
# vcf files end in: .vcf | .vcf.gz | .gvcf | .gvcf.gz | .g.vcf | .g.vcf.gz
def insert_vcf_extension(vcf_filename, extra_extension):
    VCF_EXTENSIONS = (".vcf", ".vcf.gz",
                      ".gvcf", ".gvcf.gz",
                      ".g.vcf", ".g.vcf.gz")

    extra_extension = extra_extension.strip(".")
    new_filename = None
    for ex in VCF_EXTENSIONS:
        if vcf_filename.endswith(ex):
            new_filename = vcf_filename.rstrip(ex) + \
                           "." + extra_extension + ex

    if new_filename:
        return new_filename
    else:
        return vcf_filename + extra_extension


# static method for reading sample ids from a vcf file
def get_samples_from_vcf(vcf_filename):
    vcf = VCF(vcf_filename)
    sample_ids = vcf.samples
    return sample_ids


class VCFWrapper:

    def __init__(self, *args, **kwargs):
        pass

    @classmethod
    def wrapper_factory(cls,
                        files,
                        cohort_filename,
                        use_database=False,
                        intersections=[],
                        filters=[],
                        extra_threads=0,
                        regions=[],
                        info_rules=[]):
        """
        Returns either a Database wrapper (for Wormtable)
        or a VCF parser wrapper (for CyVCF2)
        """

        if use_database:
            pass

        else:
            return CyVCFWrapper(files=files,
                                merged_filename=cohort_filename)


class CyVCFWrapper:

    def __init__(self,
                 files,
                 merged_filename,
                 filters=[],
                 extra_threads=0,
                 regions=[],
                 info_rules=[]):

        self.filenames = files
        self.cohort_filename = os.path.abspath(merged_filename)
        self.vcf = None
        self.modified_vcf = None

        # load single vcf file, assume it contains multiple samples
        if len(self.filenames) == 1:
            if not self.filenames[0] == self.cohort_filename:
                shutil.copyfile(self.filenames[0], self.cohort_filename)

        # merge multiple input vcf files
        else:
            self._merge(self.filenames,
                        self.cohort_filename,
                        threads=str(extra_threads),
                        regions=regions,
                        info_rules=info_rules)

        self.vcf = VCF(self.cohort_filename)

        # VEP Annotation parser
        self.vep_parser = VEPAnnotation(self.get_headers())
        if not self.vep_parser.vep_annotated:
            self.vep_parser = None

    def _not_indexed(self,
                     filenames):
        """
        Check if for all vcf files in filenames, an index file exists
        """
        not_indexed = []

        for file in filenames:
            index_file = file + ".tbi"

            if not os.path.exists(index_file):
                not_indexed.append(file)

            # if the index file is older than the vcf, it is removed
            # because it has to be considered outdated.
            # the vcf is then reindexed
            elif not os.path.getmtime(index_file) == os.path.getmtime(file):
                os.remove(index_file)
                not_indexed.append(file)

        return not_indexed

    def _index(self,
               filenames,
               tabix_format=True,
               reindex=False,
               threads=None):
        """
        Create VCF index files if they don't already exist by calling
        bcftools
        """

        # TODO: other index formats
        index_missing = self._not_indexed(filenames)

        if index_missing or reindex:

            # TODO: threading
            for filename in index_missing:
                command = ["bcftools", "index"]
                # input files
                command += [filename]
                # index format
                if tabix_format:
                    command += ["--tbi"]
                if threads:
                    command += ["--threads", threads]

                _ = subprocess.run(command, check=True)

    def _merge(self,
               filenames,
               outfile,
               threads=0,
               regions=None,
               info_rules=None):
        """
        Use 'vcftools/bcftools merge' to merge multiple vcf files for
        the same cohort into a single vcf.
        If vcfs are not indexed, it will throw an error.
        If bcftools is not installed, it will throw an error
        """
        check_dependency("bcftools")

        self._index(filenames)

        command = ["bcftools", "merge"]

        # input files
        command += filenames

        # extra options for bcftools
        # output file
        command += ["-o", outfile]

        # output type should be compressed vcf: .vcf.gzip
        command += ["-O", "z"]

        # Genotypes at missing types should be considerd 0/0 (Reference)
        command += ["-0"]

        # info field rules
        # 1. merge the VEP annotations by joining
        command += ["-i", "CSQ:join"]
        # 2. apply further joining operations
        # TODO

        # force rename of sampels with the same name
        command += ["--force-samples"]

        # thread control
        if threads:
            command += ["--threads", threads]
        # region filter
        if regions:
            command += ["-r", regions]

        process = subprocess.run(command, check=True)
        self.vcf = VCF(outfile)

        # add information about the merger to the new vcf file header
        infodict = {
            "ID": "bcftools-merge",
            "Number": "1",
            "Type": "Merge",
            "Description": "Merging of cohort samples into single vcf. " +
                           "Commandline: {}".format(" ".join(process.args))
        }
        self.vcf.add_info_to_header(infodict)
        self._index([outfile])

    def _merge_with_intersect(self,
                              filenames,
                              outfile,
                              threads=None,
                              regions=None,
                              info_rules=None):
        """
        Merge vcfs and compute intersections.
        """
        pass

    def __iter__(self,
                 samples=[]):
        """
        Return a variant iterator
        """
        if samples:
            vcf = VCF(self.cohort_filename, samples=samples)
        else:
            vcf = VCF(self.cohort_filename)

        return vcf

    def __next__(self,
                 samples=[]):
        """
        Iterate over the records in the underlying vcf representation
        """
        if samples:
            vcf = VCF(self.cohort_filename, samles=samples)
        for var in vcf:
            vcf = VCF(self.cohort_filename)
            yield var

    def _variant_intersection(self,
                              record,
                              alt_calls,
                              call_rate):
        """
        Returns True if the record is an intersection with a callrate of
        at least p
        """

        if ((record.num_het + record.num_hom_alt) / (record.num_called)
            >= alt_calls) and \
           record.call_rate >= call_rate:
            return True
        else:
            return False

    def intersect(self,
                  sample_groups,
                  alt_rate,
                  call_rate):
        """
        Iterates over the cohort vcf and intersects each of the sample groups.
        Creates a new VCF File with only the records that intersect for
        any of the groups.
            sample_groups       :   List of sample groups, which are itself
                                    lists containing sample IDs
        """

        # make the entire record available in each iteration step
        vcfs = [VCF(self.cohort_filename)]

        # create iterators that only return records with samples from a group
        for sample_group in sample_groups:
            vcfs.append(VCF(self.cohort_filename, samples=sample_group))

        # create a new VCF file, where only records that intersect
        # in at least one of the sample_groups are added
        extra_ext = ".isec_temp"
        new_filename = insert_vcf_extension(self.cohort_filename, extra_ext)
        intersected_vcf = Writer(new_filename, vcfs[0])

        # add info about the performed intersection to the new vcf
        info = {
            "ID": "DiscavarIntersect",
            "Description": "Intersections of {} ".format(sample_groups) +
                           "with a call rate of at least {}".format(call_rate)
        }
        intersected_vcf.add_filter_to_header(info)

        # call_rate is the percentage of variants in a record that have to
        # have the same genotype in at least one allele.
        # it therefore can't be smaller than 1.
        call_rate = max(min(1, call_rate), 0)

        for records in zip(*vcfs):

            # records[0] is the entire record with all samples.
            # records[1:] each only contain only information samples
            # from one group
            for var in records[1:]:
                if self._variant_intersection(var, alt_rate, call_rate):
                    intersected_vcf.write_record(records[0])
                    break

        # replace the current cohort file by the one with intersection applied
        os.remove(self.cohort_filename)
        os.rename(new_filename, self.cohort_filename)
        self.vcf = intersected_vcf
        self.vcf.close()

        # reindex the cohort vcf file
        self._index([self.cohort_filename], reindex=True)

    def _variant_subtraction(self,
                             minuend,
                             subtrahend,
                             call_rate1,
                             call_rate2):
        """
        Returns true if the variants in record2 are called no more than
        call_rate2 and the variants in record1 are called with a call rate
        of at least call_rate1
        """
        if (minuend.call_rate >= call_rate1) and \
           (subtrahend.call_rate < call_rate2):
            return True
        return False

    def subtract(self,
                 minuend_samples,
                 subtrahend_samples,
                 call_rate1,
                 alt_ratio1,
                 call_rate2,
                 alt_ratio2):
        """
        Subtracts the variantsin subtrahend_samples from the variants in
        minuend_samples.
        @param minuend_samples     : list of sample ids that are to represent
                                     the minuend of this operation
        @param subtrahend_samples  : list of sample ids that are to be
                                     subtracted from the first group of samples
        @param call_rate1          : percentage of samples in group1 that have
                                     to intersect.
        @param call_rate2          : percentage of samples in group2 that can
                                     at most share the genotype of the samples
                                     in group1
        """

        # make the entire record available in each iteration step
        vcfs = [VCF(self.cohort_filename)]

        # create separate iterators for the two groups
        # subtrahend
        vcfs.append(VCF(self.cohort_filename, samples=minuend_samples))
        vcfs.append(VCF(self.cohort_filename, samples=subtrahend_samples))

        # create a new VCF file to which only records that appear in the
        # minuend group uniquely are added to
        extra_ext = ".subt"
        new_filename = insert_vcf_extension(self.cohort_filename, extra_ext)
        subtracted_vcf = Writer(new_filename, vcfs[0])
        # add info about the performed intersection to the new vcf
        info = {
            "ID": "Intersection",
            "Description": "Subtraction of the samples {} with a call rate " +
                           "of at most {} "
                           .format(subtrahend_samples, call_rate2) +
                           "from the samples {} with a call rate of " +
                           "at least {} for every record"
                           .format(minuend_samples, call_rate1)
        }
        subtracted_vcf.add_filter_to_header(info)

        call_rate1 = min(1, call_rate1)
        call_rate2 = min(1, call_rate2)

        for orig, minu, subt in zip(*vcfs):

            # orig is the entire record with both minuend and subtrahend
            # samples
            # rminu is the record with only minuend samples
            # subt is the record with only subtrahend samples
            if self._variant_subtraction(minu, subt, call_rate1, call_rate2):
                subtracted_vcf.write_record(orig)

        # replace the old cohort file
        self.cohort_filename = new_filename
        self.vcf = subtracted_vcf
        self.vcf.close()

        # reindex the cohort vcf file
        self._index([self.cohort_filename], reindex=True)

    def get_headers(self):
        """
        Return all headers of the vcf
        """
        self.vcf = VCF(self.cohort_filename)
        headers = list(self.vcf.header_iter())
        return headers

    def get_header_iter(self):
        """
        Return an iterator for all headers of the vcf
        """
        self.vcf = VCF(self.cohort_filename)
        return self.vcf.header_iter()

    def get_variants_at(self, regions=[]):
        """
        Returns all variants that lie in the given regions
        """
        self.vcf = VCF(self.cohort_filename)
        variants = []

        for region in regions:
            for var in self.vcf(region):
                variants += var

        return variants

    def get_sample_names(self):
        """
        Return the column header for all sample columns in the
        underlying vcf file
        """
        self.vcf = VCF(self.cohort_filename)
        return self.vcf.samples

    # methods for accessing certain variant attributes based
    # on variant types used by instances of this class
    @classmethod
    def get_var_chrom(cls, var):
        """
        return the chromosome for the record
        """
        return var.CHROM

    @classmethod
    def get_var_pos(cls, var):
        """
        Return the position of the record on a chromosome
        """
        return var.POS

    @classmethod
    def get_var_gq(cls, var):
        """
        return genotype qualities for all samples in a record
        """
        return var.QUAL

    @classmethod
    def get_var_ref(cls, var):
        """
        Return the reference genotype
        """
        return var.REF

    @classmethod
    def get_var_alt(cls, var):
        """
        Return alternative genotype for a var
        """
        return var.alt

    @classmethod
    def get_var_gt(cls, var):
        """
        Return bases for all samples
        """
        return var.gt_bases

    @classmethod
    def get_var_gtcode(cls, var):
        """
        Returns an array with genotype codes for each sample
        CyVCF2 Convention:
            HOM_REF=0,
            HET=1.
            For gts012=True
            HOM_ALT=2,
            UKNOWN=3
        """
        return var.gt_types

    @classmethod
    def get_var_type(cls, var):
        """
        Return the type of the variant
        """
        return var.var_type

    @classmethod
    def get_var_filter(cls, var):
        """
        Returns the value of the filter column for a variant
        Possible values are: 'PASS' or a comma separated list of filters
        the variant did NOT pass.
        """
        return var.FILTER

    @classmethod
    def get_var_info(cls, var, field=""):
        """
        Return a dictionary object
        """
        if field:
            info = var.INFO.get(field, "")
        else:
            info = [i for i in var.INFO]

        return info

    # methods for making queries
    def apply(self,
              filters=[]):
        """
        Applies a set of filters all at once to the underlying vcf,
        adds a header line and replaces it with the current vcf.
        """
        if filters:

            # reopen the vcf iterator
            self.vcf = VCF(self.cohort_filename)

            # new file where only the records that pass all filters are written
            # to
            new_filename = insert_vcf_extension(self.cohort_filename, "temp")
            new_vcf = Writer(new_filename, self.vcf)

            passed_records = 0
            for record in self.vcf:
                all_filters_passed = True

                for fltr in filters:
                    if not fltr.passes(record, self):
                        all_filters_passed = False
                        break

                if all_filters_passed:
                    passed_records += 1
                    new_vcf.write_record(record)

            print("Records passed: {}".format(passed_records))

            # add filter information to the vcf header
            for fltr in filters:
                new_vcf.add_filter_to_header(fltr.get_filter_info())

            # replace the current cohort file by the one with filters applied
            os.remove(self.cohort_filename)
            os.rename(new_filename, self.cohort_filename)
            self.vcf = new_vcf
            self.vcf.close()

            # reindex the cohort vcf file
            self._index([self.cohort_filename], reindex=True)

    def query(self,
              filters=[]):
        """
        Queries the underlying vcf for variants that satisfy given
        intersection and filter rules as well as direct queries.
        Returns an iterator objects that allows loping over the variant
        records that result from the query
        """
        if filters:
            self.vcf = VCF(self.cohort_filename)

            for record in self.vcf:
                all_filters_passed = True

                for fltr in filters:
                    if not fltr.passes(record, self):
                        all_filters_passed = False
                        break

                if all_filters_passed:
                    yield record

    def to_vcf(self):
        """
        Save the underlying vcf representation
        """
        self.vcf.close()


# -------------- #
#  FUTURE WORK   #
# -------------- #
class WormtableWrapper:
    """
    NOT IMPLEMENTED
    """
    def __init__(self, *args, **kwargs):
        self.wormtable = None

    def get_headers(self):
        """
        Return all headers of the vcf
        """
        pass

    def get_header_iter(self):
        """
        Return an iterator for all headers of the vcf
        """
        pass

    def get_variants_at(self, regions=[]):
        """
        Returns all variants that lie in the given regions
        """
        pass

    def apply(self,
              intersection=None,
              query=None,
              filter=None,
              intersections=[],
              queries=[],
              filters=[]):
        """
        Applies a set of intersection rules, queries and filters all
        at once to the underlying wormtable, adds a header line and replaces
        it with the current vcf.
        """
        pass

    def query(self,
              intersection=None,
              query=None,
              filter=None,
              intersections=[],
              queries=[],
              filters=[]):
        """
        Queries the underlying wormtable for variants that satisfy given
        intersection and filter rules as well as direct queries
        """
        pass

    def to_vcf(self):
        """
        Save the underlying wormtable to a vcf
        """
        pass
