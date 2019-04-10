from .cli import warning
import re


class VEPAnnotation:
    """
    Takes a vcf header and looks for VEP annotations. If the vcf file was
    annotated with VEP, an INFO header with id "CSQ" is assumed to exist.


    Currently supports VEP Version v93
    VEP Anotations in INFO column of VCF files:
        Consequence annotations predict the effect of variants.
        ID: CSQ
        Format:
            Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|
            Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|
            CDS_position|Protein_position|Amino_acids|Codons|
            Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|
            SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|APPRIS|CCDS|ENSP|
            SWISSPROT|TREMBL|UNIPARC|GENE_PHENO|SIFT|PolyPhen|
            DOMAINS|miRNA|AF|AFR_AF|AMR_AF|EAS_AF|EUR_AF|SAS_AF|
            AA_AF|EA_AF|gnomAD_AF|gnomAD_AFR_AF|gnomAD_AMR_AF|
            gnomAD_ASJ_AF|gnomAD_EAS_AF|gnomAD_FIN_AF|gnomAD_NFE_AF|
            gnomAD_OTH_AF|gnomAD_SAS_AF|MAX_AF|MAX_AF_POPS|FREQS|
            CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|
            HIGH_INF_POS|MOTIF_SCORE_CHANGE
    """

    # separator for fields in the INFO column of a vcf file
    COLUMN_SEPARATOR = ","

    # the separator for fields of a VEP Annotation
    SEPARATOR = "|"
    SUB_SEPARATOR = "&"
    INFO_COL_NAME = "CSQ"

    def __init__(self, vcf_headers):
        self.vep_annotated = False
        self.fields = self.get_ann_fields(vcf_headers)

        if not self.vep_annotated:
            warning("VCF not VEP Annotated!")

    def get_ann_fields(self, headers):
        """
        Returns a list of fields that are to be expected in the CSQ
        INFO column annotated by VEP.
        VEP Annotation is determined by checking for an INFO header
        with ID=CSQ.
        """

        for header in headers:

            # cyvcf header representation
            if header.type == "INFO":

                # get a python dict with all info fields in the header line
                info = header.info()
                info_id = info.get("ID", "")
                description = info.get("Description", "")

                if description and info_id == "CSQ":

                    # find startposition of the fields description
                    fields_pattern = re.compile(
                        r'(\w+' + re.escape(self.SEPARATOR) + r')+\w+')
                    match = fields_pattern.search(description)
                    if match:
                        start = match.start()
                        fields = description[start:].split(self.SEPARATOR)
                        self.vep_annotated = True

                        return fields
        return None

    def parse_annotation(self, ann_string, separator=""):
        """
        Creates a dictionary with(fieldname, variant_annotation)
        pairs for the VEP annotation of a every variant in a given record
        and returns them as a list.

        @param ann_string   : A string representing the VEP Annotation Column
                              of a record in a VCF file
        @param separator    : A string representing the character used to
                              separate annotation fields for a variant.
        """
        if self.vep_annotated:

            # check if a different seperator should be used
            if not separator:
                separator = self.SEPARATOR

            # fetch the annotations for every separate variant in a record
            record_annotations = ann_string.split(self.COLUMN_SEPARATOR)

            # list to save the dict with <field>: <value> pairs for
            # each variant
            var_annotations = []

            for var_annotation_string in record_annotations:

                # obtain the values for each VEP Annotation Field 
                # (see class comment)
                var_annotation_values = ann_string.split(separator)

                # construct a dictionary object for each variant annotation
                var_annotation = dict()

                for k, ann in zip(self.fields, var_annotation_values):

                    # check if the annotation value is a list of values
                    if ann.find(self.SUB_SEPARATOR) != -1:
                        var_annotation[k] = ann.split(self.SUB_SEPARATOR)
                    else:
                        var_annotation[k] = ann

                var_annotations.append(var_annotation)

            return var_annotations

        return None
