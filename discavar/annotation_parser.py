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
    # the separator for fields of the annotation string
    SEPARATOR = "|"
    SUB_SEPARATOR = "&"
    INFO_COL_NAME = "CSQ"

    def __init__(self, vcf_headers):
        self.vep_annotated = False
        self.fields = self.get_ann_fields(vcf_headers)

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
        Returns a dictionary of (fieldname, variant_annotation) pairs for
        the VEP annotation of a variant
        """
        if not separator:
            separator = self.SEPARATOR

        values = ann_string.split(separator)
        annotations = dict(zip(self.fields, values))
        #for k, ann in annotations:
        #    annotations[k] = ann.split(self.SUB_SEPARATOR)
        return annotations
