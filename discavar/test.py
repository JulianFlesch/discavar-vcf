from discavar.analysis import Cohort
from discavar.annotation_parser import VEPAnnotation
from discavar.filters import AnnotationFilter, VariantFilter
from discavar.vcfwrapper import CyVCFWrapper
import os
import time
tsvfile = "/home/julian/Uni/Bachelorarbeit/vcfs/input.tsv"
os.chdir("../../vcfs")

# t1 = time.time()

cohorts = Cohort.from_tsv(tsvfile, build_cohort=False)
c = cohorts[0]
c.variants = CyVCFWrapper(["mrkh_analysis_vcf2/cohort_1.vcf.gz"], "")
#c.build()

var_filter = VariantFilter(min_gq=20)

ann_parser = VEPAnnotation(c.variants.get_headers())

consequences = ["stop_lost", "stop_gained", "start_gained", "start_lost",
                "frameshift_variant", "missense_variant",
                "NMD_transcript_variant", "coding_sequence_variant"]

impacts = ["MODERATE", "HIGH", "MODIFIER"]

# significance = ["unknown_significance", "pathogenic", "risk_factor", "association"]
# significance_not = ["benign", "likely_benign"]
ann_filter = AnnotationFilter(ann_parser,
#                             Consequence__in=consequences,
                              IMPACT__in=impacts,
#                              CLIN_SIG__in=significance
                              )

#c.add_filter(var_filter)
#c.add_filter(ann_filter)
#c.intersect(0.5, 0.5)
# t3 = time.time()
#c.apply_filters()
# t4 = time.time()

# _ = c.manhattan()

# _ = c.vtypes_dist()
