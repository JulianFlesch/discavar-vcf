#!/home/julian/Uni/Bachelorarbeit/ba-venv/bin/python3

"""
    DISCAVAR: VCF FILE ANALYTICS
    ============================

    An variation analysis tool developed as part of a bachelor thesis at the 
    University of Tuebingen.

    Contacts:

        Author:         Julian Flesch
                        julian.flesch@student.uni-tuebingen.de

        Supervisor:     Nico Weber
        Project Owner:  Prof. Dr. Oliver Kohlbacher

    Thesis:             Discavar-VCF detecting disease causing variant
                        in familial NGS Data with Nextflow pipelines.

    Documentation:      TODO

    Licence:            TODO
"""

from cli import parse_args
from analysis import AnalysisControl


def main():
    # parse the command line arguments and save it in a namespace
    job = parse_args()

    # use the arguments defined in job for the right subcommand
    cmd = job.subcommand

    analysis_control = AnalysisControl(analysis=cmd,
                                       job=job)

    analysis_control.run()


# for testing
if __name__ == "__main__":
    main()

    # TODO: deployed version should be called via a startscript in bin/
