from cyvcf2 import VCF
from cli import warning, error, Description


# static method for applying a check multiple samples
def all_pass(variants, check):
    for var in variants:
        if not check(var):
            return False
    return True


class AnnotationFilter:
    """
    AnnotationFilter takes an annotation parser and a list of keyword arguments
    in the form FIELDNAME__op where "__op" is one of:
        __gt        :   greater than a given number,
        __ge        :   greater or equal than a given number,
        __lt        :   less than a given number,
        __le        :   less or equal to a given number,
        __eq        :   equal to a string,
        __in        :   in a list list of strings, or
        __notin     :   not in a list of strings
        in the annotation FIELDNAME
    """
    def __init__(self,
                 annotation_parser,
                 **kwargs):

        # keyword args have to end in exactly one of these:
        self.gt_suffix = "__gt"
        self.ge_suffix = "__ge"
        self.lt_suffix = "__lt"
        self.le_suffix = "__le"
        self.in_suffix = "__in"
        self.notin_suffix = "__notin"
        self.eq_suffix = "__eq"

        self.annotation_parser = annotation_parser
        self.gt_rules = []
        self.ge_rules = []
        self.lt_rules = []
        self.le_rules = []
        self.eq_rules = []
        self.in_rules = []
        self.notin_rules = []

        for key, value in kwargs.items():
            if key.split("__")[0] not in self.annotation_parser.fields:
                warning("Unknown Annotation Field")

            if key.endswith(self.gt_suffix):
                self.gt_rules.append((key.rstrip(self.gt_suffix), value))
            elif key.endswith(self.ge_suffix):
                self.ge_rules.append((key.rstrip(self.ge_suffix), value))
            elif key.endswith(self.lt_suffix):
                self.lt_rules.append((key.rstrip(self.lt_suffix), value))
            elif key.endswith(self.le_suffix):
                self.le_rules.append((key.rstrip(self.le_suffix), value))
            elif key.endswith(self.eq_suffix):
                self.eq_rules.append((key.rstrip(self.eq_suffix), value))
            elif key.endswith(self.in_suffix):
                self.in_rules.append((key.rstrip(self.in_suffix), value))
            elif key.endswith(self.notin_suffix):
                self.notin_rules.append((key.rstrip(self.notin_suffix), value))
            else:
                warning("Skipping: Invalid Filter operation.")

    def get_filter_info(self):
        """
        Returns a dict with the keys ID and Description.
        """
        filter_id = "AnnKW"
        filter_description = "Annotation Keyword Filter {} Version {}".format(
                    Description.NAME, Description.VERSION)

        return {"ID": filter_id, "Description": filter_description}

    def passes(self, record, wrapper_cls):
        # get the annotation string from the record returned by a vcf_wrapper
        ann_string = wrapper_cls.get_var_info(
                        record,
                        field=self.annotation_parser.INFO_COL_NAME)
        # parse the annotation to obtain a dictionary object
        annotations = self.annotation_parser.parse_annotation(ann_string)

        # check for the satisfaction of all specified rules
        # greater than(gt), less than(lt), ... rules are for numeric operations
        # the annotation field must therefor be a number
        for field, gt_value in self.gt_rules:
            if not int(annotations.get(field, 0)) > gt_value:
                return False

        for field, ge_value in self.ge_rules:
            if not int(annotations.get(field, 0)) >= ge_value:
                return False

        for field, lt_value in self.lt_rules:
            if not int(annotations.get(field, lt_value)) < lt_value:
                return False

        for field, le_value in self.le_rules:
            if not int(annotations.get(field, lt_value+1)) <= le_value:
                return False

        # For "eq", "in" and "notin" rules the annotation field is expected to
        # potentially have more than one value(="SUB_FIELD").
        # The annotation parser has a subfield-separator value, by which the
        # annotation can be further divided.
        for field, eq_value in self.eq_rules:
            for ann in annotations.get(field, "").split(
                 self.annotation_parser.SUB_SEPARATOR):
                if ann == eq_value:
                    break
            else:
                return False

        for field, in_values in self.in_rules:
            for ann in annotations.get(field, "").split(
                 self.annotation_parser.SUB_SEPARATOR):
                if ann in in_values:
                    break
            else:
                return False

        for field, notin_values in self.notin_rules:
            for ann in annotations.get(field, "").split(
                 self.annotation_parser.SUB_SEPARATOR):
                if ann in notin_values:
                    return False

        # Variant has passed all filters!
        return True


class VariantFilter:
    """
    VariantFilter is a filter for variants based on a set of parameters.
    These are:
        min_gq      :   genotype quality
        genes       :   A set of genes that should exclusively be
                        considered. Note: Only works with annotations!
        excl_genes  :   A set of genes that should be excluded.
                        Note: Only works with annotations!
        intervals   :   A set of intervals that sould exclusively be
                        considered
        excl_intervals: A set of itervals that sould be excluded
        vtypes      :   A set of variation types that should exclusively
                        be considered.
        excl_vtypes :   A set of variation types that should be ignored
        perc_pass   :   A float between 0 and 1 describing the percentage
                        of genotype calls in multi-sample vcfs that have to
                        pass the filter.
    """
    def __init__(self,
                 min_gq=0,
                 min_mq=0,
                 min_dp=1,
                 genes=[],
                 excl_genes=[],
                 intervals=[],
                 excl_intervals=[],
                 vtypes=[],
                 excl_vtypes=[],
                 perc_pass=0.5):

        if min_gq:
            try:
                self.min_gq = min_gq[0]
            except TypeError:
                try:
                    self.min_gq = float(min_gq)
                except ValueError:
                    warning("Invalid min-gq value. Using default ..")
                    self.min_gq = 0
        else:
            self.min_gq = 0

        if min_mq:
            try:
                self.min_mq = min_mq[0]
            except TypeError:
                try:
                    self.min_mq = float(min_mq)
                except ValueError:
                    warning("Invalid min-mq value. Using default ..")
                    self.min_mq = 0
        else:
            self.min_mq = 0

        if min_dp:
            try:
                self.min_dp = min_dp[0]
            except TypeError:
                try:
                    self.min_dp = float(min_dp)
                except ValueError:
                    warning("Invalid min-qc value. Using default ..")
                    self.min_dp = 0
        else:
            self.min_dp = 1

        if perc_pass:
            self.perc_pass = perc_pass
        else:
            perc_pass = 1

        self.genes = genes
        self.excl_genes = excl_genes
        self.intervals = intervals
        self.excl_intervals = excl_intervals
        self.vtypes = vtypes
        self.excl_vtypes = excl_vtypes

    def __str__(self):
        """
        Returns a descriptive string of the filter.
        Note: Output should not be considered stable (e.g. for use in
        VCF headers)
        """
        header = '<VariantFilter, {}>'
        header_id, header_descr = self.get_filter_info()
        return header.format(header_descr)

    def get_filter_info(self):
        """
        Returns a dict with the keys ID and Description.
        """
        filter_id = "{} Version={}".format(
                    Description.NAME, Description.VERSION)
        filter_description = "Genotype Quality, Variant Type, " + \
                             "Chromosome Region, and Gene Filter"

        return {"ID": filter_id, "Description": filter_description}

    def passes(self, variant, wrapper_cls, annotations=None):

        # genotype quality filter
        if wrapper_cls.get_var_gq(variant) < self.min_gq:
            return False

        # read depth filter
        if wrapper_cls.get_var_info(variant, field="DP") < self.min_dp:
            return False

        # read mapping quality
        if wrapper_cls.get_var_info(variant, field="MQ") < self.min_mq:
            return False

        # filter by variant type
        if wrapper_cls.get_var_type(variant) in self.excl_vtypes:
            return False

        if wrapper_cls.get_var_type(variant) not in self.vtypes:
            return False

        # filter by chromosomal region
        for excl in self.excl_intervals:
            if wrapper_cls.get_var_pos(variant) > excl.start and \
               wrapper_cls.get_var_pos(variant) < excl.end and \
               wrapper_cls.get_var_chrom(variant) == excl.chrom:
                return False

        for incl in self.intervals:
            if wrapper_cls.get_var_pos(variant) > incl.start and \
               wrapper_cls.get_var_pos(variant) < incl.end and \
               wrapper_cls.get_var_chrom(variant) == incl.chrom:
                break
        else:
            return False

        # all filters passed
        return True
