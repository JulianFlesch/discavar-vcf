from analysis import Cohort
from annotation_parser import VEPAnnotation
from filters import AnnotationFilter, VariantFilter
from vcfwrapper import CyVCFWrapper
import os
import time
tsvfile = "/home/julian/Uni/Bachelorarbeit/vcfs/test_input.tsv"
os.chdir("../../vcfs")

# t1 = time.time()

cohorts = Cohort.from_tsv(tsvfile, build_cohort=False)
c = cohorts[0]
#c.variants = CyVCFWrapper(["cohort_1.vcf.gz"], "")
c.build()

var_filter = VariantFilter(min_gq=20)

ann_parser = VEPAnnotation(c.variants.get_headers())

# consequences = ["stop_lost", "stop_gained", "start_gained", "start_lost",
#                "frameshift_variant", "missense_variant",
#                "NMD_transcript_variant", "coding_sequence_variant"]

impacts = ["MODERATE", "HIGH", "MODIFIER"]

# significance = ["benign", "likely_benign"]

ann_filter = AnnotationFilter(ann_parser,
#                             Consequence__in=consequences,
                              IMPACT__in=impacts,
#                              CLIN_SIG__in=significance
                              )

c.add_filter(var_filter)
c.add_filter(ann_filter)
c.intersect(0.5, 0.5)
# t3 = time.time()
c.apply_filters()
# t4 = time.time()

# _ = c.manhattan()

# _ = c.vtypes_dist()
