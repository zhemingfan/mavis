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


def setup_library_folders(library):
    """
    Set up the initial directories
    """
    def mkdir_relative(*args):
        return _util.mkdirp(os.path.join(config['output_dir'], *args))

    mkdir_relative(library, 'cluster')
    mkdir_relative(library, 'validate')
    mkdir_relative(library, 'annotate')
    stamp_file = os.path.join(config['output_dir'], library, 'SETUP.COMPLETE')
    with open(stamp_file, 'w') as fh:
        fh.write('')


def setup_pipeline_folder():
    """
    Set up the initial directories
    """
    def mkdir_relative(*args):
        return _util.mkdirp(os.path.join(config['output_dir'], *args))
    validate_config(config, bam_stats=True)
    mkdir_relative('converted_outputs')
    mkdir_relative('pairing')
    mkdir_relative('summary')

    config_path = os.path.join(config['output_dir'], 'mavis.config.json')
    _util.LOG(f'writing: {config_path}')
    with open(config_path, 'w') as fh:
        fh.write(json.dumps(config, sort_keys=True, indent=4))

    for library in sorted(list(config['libraries'])):
        setup_library_folders(library)
    stamp_file = os.path.join(config['output_dir'],'SETUP.COMPLETE')
    with open(stamp_file, 'w') as fh:
        fh.write('')


def skip_validate(config):
    return 'validate' in config.get('skip_stage', [])

mavis_exe = shutil.which('mavis')

def output_dir(*paths):
    return os.path.join(config['output_dir'], *paths)


libraries = sorted(list(config['libraries']))

config_path = output_dir('mavis.config.json')
if not os.path.exists(config_path):
    setup_pipeline_folder()

ruleorder: convert > cluster

# "all" rule MUST be the first rule
rule all:
    input: expand(output_dir('{library}/cluster/MAVIS.COMPLETE'), library=libraries)


rule convert:
    output: output_dir('converted_outputs/{alias}.tab')
    params:
        file_type=lambda w: config['convert'][w.alias]['file_type'],
        strand_specific=lambda w: config['convert'][w.alias]['strand_specific'],
        assume_no_untemplated=lambda w: config['convert'][w.alias]['assume_no_untemplated'],
        input_files=lambda w: config['convert'][w.alias]['inputs']
    shell:
        'mavis convert --file_type {params.file_type}'
            + ' --strand_specific {params.strand_specific}'
            + ' --assume_no_untemplated {params.assume_no_untemplated}'
            + ' --inputs {params.input_files}'
            + ' --outputfile {output}'
            + ' &> ' + output_dir('converted_outputs/snakemake.{wildcards.alias}.log.txt; echo $?')


def get_cluster_inputs(w):
    print('cluster inputs', w.library, config['libraries'][w.library]['assign'], file=sys.stderr)
    return  config['libraries'][w.library]['assign']


# def get_cluster_inputs(wildcards):
#     all_files = []
#     for library in libraries:
#         for checkpoint_output in checkpoints.library_setup.get(library=library, **wildcards).output:
#             print('checkpoint_output', checkpoint_output)
#             target_path = os.path.join(checkpoint_output, 'batch-{job_id}.tab')
#             for job_id in glob_wildcards(target_path).job_id:
#                 all_files.append(os.path.join(checkpoint_output.replace('/cluster', '/annotate'), f'batch-{job_id}', 'annotations.tab'))
#     print('aggregate_jobs', all_files)
#     return all_files

rule cluster:
    input: get_cluster_inputs
    output: output_dir('{library}/cluster/MAVIS.COMPLETE')
    shell:
        'mavis cluster --config ' + config_path
            + ' --library {wildcards.library}'
            + ' --inputs {input}'
            + ' --output ' + output_dir('{wildcards.library}/cluster')
            + ' &> ' + output_dir('{wildcards.library}/snakemake.cluster.log.txt; echo $?')


# rule validate:
#     input: expand(output_dir('{library}/cluster/batch-{{job_id}}.tab'), library=libraries)
#     params:
#         dirname=output_dir('{library}/validate/batch-{job_id}')
#     output: stamp=output_dir('{library}/validate/batch-{job_id}/MAVIS.COMPLETE'),
#         result=output_dir('{library}/validate/batch-{job_id}/validation-passed.tab')
#     shell:
#         'mavis validate --config ' + config_path
#             + ' --library {wildcards.library}'
#             + ' --inputs {input}'
#             + ' --output {params.dirname}'
#             + ' &> ' + output_dir('{wildcards.library}/validate/snakemake.batch-{wildcards.job_id}.log.txt')


# rule annotate:
#     input: expand('{library}/validate/batch-{{job_id}}/validation-passed.tab', library=libraries)
#     output: stamp=output_dir('{library}/annotate/batch-{job_id}/MAVIS.COMPLETE'),
#         result=output_dir('{library}/annotate/batch-{job_id}/annotations.tab')
#     shell:
#         'mavis annotate --config ' + config_path
#             + ' --library {wildcards.library}'
#             + ' --inputs {input}'
#             + ' --output {wildcards.library}/annotate/batch-{wildcards.job_id}'
#             + ' &> ' + output_dir('{wildcards.library}/annotate/snakemake.batch-{wildcards.job_id}.log.txt')


# def aggregate_jobs(wildcards):
#     '''
#     aggregate the file names of the random number of files
#     generated at the cluster step
#     '''
#     w = {'library': libraries}
#     w.update(wildcards)
#     all_files = []
#     # patts = expand(output_dir('{library}/annotate', 'batch-{{job_id}}.tab'), library=libraries)
#     # print(patts)
#     # for patt in patts:
#     #     globbed = glob_wildcards(patt)
#     #     all_files.extend(expand(patt, job_id=globbed.job_id))
#     # print(globbed)
#     for lib in libraries:
#         for checkpoint_output in checkpoints.cluster.get(library=lib, **wildcards).output:
#             print('checkpoint_output', checkpoint_output)
#             target_path = os.path.join(checkpoint_output, 'batch-{job_id}.tab')
#             for job_id in glob_wildcards(target_path).job_id:
#                 all_files.append(os.path.join(checkpoint_output.replace('/cluster', '/annotate'), f'batch-{job_id}', 'annotations.tab'))
#     print('aggregate_jobs', all_files)
#     return all_files


# rule pairing:
#     input: aggregate_jobs
#     output: stamp=output_dir('pairing/MAVIS.COMPLETE'),
#         result=output_dir('pairing/mavis_paired_' + '_'.join(libraries) + '.tab')
#     params:
#         dirname=output_dir('pairing')
#     shell:
#         'mavis pairing --config ' + config_path
#             + ' --inputs {input}'
#             + ' --output {params.dirname}'
#             + ' &> ' + output_dir('snakemake.pairing.log.txt')

# rule summary:
#     input: rules.pairing.output.result,
#     output: stamp=output_dir('summary/MAVIS.COMPLETE')
#     params:
#         dirname=output_dir('summary')
#     shell:
#         'mavis summary --config ' + config_path
#             + ' --inputs {input}'
#             + ' --output {params.dirname}'
#             + ' &> ' + output_dir('snakemake.summary.log.txt')
