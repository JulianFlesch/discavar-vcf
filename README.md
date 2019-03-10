# DISCAVAR VCF
VCF utility written in Python for calldata and annotation based filtering, interactive exploration and identification of disease causative variations.

# Installation

## GitHub

`pip install -r requirements.txt`

## PIP

`pip install discavar`

## Anaconda

TODO

## Dockerhub

TODO

# Dependencies

For merging different vcf files into single cohort file, `bcftools` is required.

# Running

From commandline with the executable `discavar`. See `discavar --help for details`

Interactively form a python shell by importing from discavar. See documentation for details.

e.g.: 
``` from discavar import Cohort
    
    sample_file = "cohorts_input.tsv"
    cohorts = Cohort.from_tsv(sample_file, build_cohorts=True)
    
    # subtraction analysis
    for cohort in cohorts:
        cohort.healthy_vs_diseased()

    # calldata filter
    var_filter = VariantFilter(
                            min_gq=20,
                            min_dp=10,
                            excl_intervals=["chr1"])

    # annotation filter
    ann_parser = VEPAnnotation(c.variants.get_headers())
    impacts = ["MODERATE", "HIGH", "MODIFIER"]
    significance = ["benign", "likely_benign"]
    ann_filter = AnnotationFilter(ann_parser,
                                  IMPACT__in=impacts,
                                  CLIN_SIG__in=significance)

    for cohort in cohorts:
        cohort.add_filters([var_filter, ann_filter])
        cohort.apply_filters()


    # plotting
    for cohort in cohorts:
        cohort.manhattan()
        cohort.vtypes_dist()

    # cluster analysis
    ...
    # custom queries
    ...
```


