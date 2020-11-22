from snakemake.utils import validate
import os
from typing import List
from glob import glob
import shutil
import pandas
import re

from mavis.config import validate_config
from mavis import util as _util
from mavis.constants import SUBCOMMAND

validate_config(config)

def output_dir(*paths):
    return os.path.join(config['output_dir'], *paths)


libraries = sorted(list(config['libraries']))
config_path = output_dir('mavis.config.json')

jobs = glob_wildcards(output_dir('{library}/cluster/batch-{job_id}.tab'))

# "all" rule MUST be the first rule
rule all:
    input: output_dir('summary/MAVIS.COMPLETE')

rule validate:
    input: output_dir('{library}/cluster/batch-{job_id}.tab')
    params:
        dirname=lambda w: output_dir(f'{w.library}/validate/batch-{w.job_id}')
    output: output_dir('{library}/validate/batch-{job_id}/validation-passed.tab')
    shell:
        'mavis validate --config ' + config_path
            + ' --library {wildcards.library}'
            + ' --inputs {input}'
            + ' --output {params.dirname}'
            + ' &> ' + output_dir('{wildcards.library}/validate/snakemake.batch-{wildcards.job_id}.log.txt')

rule annotate:
    input: output_dir('{library}/validate/batch-{job_id}/validation-passed.tab')
    output: stamp=output_dir('{library}/annotate/batch-{job_id}/MAVIS.COMPLETE'),
        result=output_dir('{library}/annotate/batch-{job_id}/annotations.tab')
    shell:
        'mavis annotate --config ' + config_path
            + ' --library {wildcards.library}'
            + ' --inputs {input}'
            + ' --output ' + output_dir('{wildcards.library}/annotate/batch-{wildcards.job_id}')
            + ' &> ' + output_dir('{wildcards.library}/annotate/snakemake.batch-{wildcards.job_id}.log.txt')

print('expected pairing inputs', expand(rules.annotate.output.result, zip, library=jobs.library, job_id=jobs.job_id))
rule pairing:
    input: expand(rules.annotate.output.result, zip, library=jobs.library, job_id=jobs.job_id)
    output: stamp=output_dir('pairing/MAVIS.COMPLETE'),
        result=output_dir('pairing/mavis_paired.tab')
    params:
        dirname=output_dir('pairing')
    shell:
        'mavis pairing --config ' + config_path
            + ' --inputs {input}'
            + ' --output {params.dirname}'
            + ' &> ' + output_dir('snakemake.pairing.log.txt')

rule summary:
    input: rules.pairing.output.result,
    output: stamp=output_dir('summary/MAVIS.COMPLETE')
    params:
        dirname=output_dir('summary')
    shell:
        'mavis summary --config ' + config_path
            + ' --inputs {input}'
            + ' --output {params.dirname}'
            + ' &> ' + output_dir('snakemake.summary.log.txt')
