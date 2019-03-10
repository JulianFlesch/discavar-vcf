from filters import AnnotationFilter
from vcfwrapper import CyVCFWrapper as cls
from cyvcf2 import VCF
from annotation_parser import VEPAnnotation
import os

os.chdir("../../vcfs")
vcf = VCF("cohort_1.vcf.gz")
parser = VEPAnnotation(list(vcf.header_iter()))

consequences = ["stop_lost", "stop_gained", "start_gained", "start_lost",
                "3_prime_UTR_variant", "splice_donor_variant",
                "frameshift_variant", "splice_region_variant",
                "missense_variant", "NMD_transcript_variant",
                "TFBS_ablation", "TF_binding_site_variant",
                "coding_sequence_variant"]

fltr = AnnotationFilter(parser, Consequence__in=consequences)

for var in vcf:
    if fltr.passes(var, cls):
        print("one passed!!!")
        break
