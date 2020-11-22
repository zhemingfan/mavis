#!python
import argparse
import json
import logging
import os
import platform
import sys
import time
from typing import Dict

import tab

from . import __version__
from . import annotate as _annotate
from . import config as _config
from . import util as _util
from .align import get_aligner_version
from .annotate import main as annotate_main
from .cluster import main as cluster_main
from .constants import SUBCOMMAND
from .error import DrawingFitError
from .illustrate.constants import DiagramSettings
from .illustrate.diagram import draw_multi_transcript_overlay
from .illustrate.scatter import bam_to_scatter
from .pairing import main as pairing_main
from .summary import main as summary_main
from .tools import SUPPORTED_TOOL, convert_tool_output
from .util import filepath
from .validate import main as validate_main


def check_overlay_args(args, parser):
    """
    parse the overlay options and check the formatting
    """
    # check complex options
    for marker in args.markers:
        if len(marker) < 3:
            marker.append(marker[-1])
        try:
            marker[1] = int(marker[1])
            marker[2] = int(marker[2])
        except ValueError:
            parser.error('argument --marker: start and end must be integers: {}'.format(marker))

    defaults = [None, None, 0.5, None, True]
    bam_file, density, ymax, stranded = range(1, 5)

    for plot in args.read_depth_plots:
        for i, d in enumerate(defaults):
            if i >= len(plot):
                plot.append(d)
        if not os.path.exists(plot[bam_file]):
            parser.error(
                'argument --read_depth_plots: the bam file given does not exist: {}'.format(
                    plot[bam_file]
                )
            )
        try:
            plot[density] = float(plot[density])
            if plot[density] < 0 or plot[density] > 1:
                raise ValueError()
        except ValueError:
            parser.error(
                'argument --read_depth_plots: density must be an float between 0 and 1: {}'.format(
                    plot[density]
                )
            )
        try:
            if str(plot[ymax]).lower() in ['null', 'none']:
                plot[ymax] = None
            else:
                plot[ymax] = int(plot[ymax])
        except ValueError:
            parser.error(
                'argument --read_depth_plots: ymax must be an integer: {}'.format(plot[ymax])
            )
        try:
            plot[stranded] = tab.cast_boolean(plot[stranded])
        except TypeError:
            parser.error(
                'argument --read_depth_plots: stranded must be an boolean: {}'.format(
                    plot[stranded]
                )
            )
    return args


def overlay_main(
    gene_name,
    output,
    buffer_length,
    read_depth_plots,
    markers,
    annotations,
    drawing_width_iter_increase,
    max_drawing_retries,
    min_mapping_quality,
    ymax_color='#FF0000',
    **kwargs,
):
    """
    generates an overlay diagram
    """
    annotations.load()
    # check options formatting
    gene_to_draw = None

    for chrom in annotations.content:
        for gene in annotations.content[chrom]:
            if gene_name in gene.aliases or gene_name == gene.name:
                gene_to_draw = gene
                _util.LOG(
                    'Found target gene: {}(aka. {}) {}:{}-{}'.format(
                        gene.name, gene.aliases, gene.chr, gene.start, gene.end
                    )
                )
                break
    if gene_to_draw is None:
        raise KeyError('Could not find gene alias or id in annotations file', gene_name)

    settings = DiagramSettings(**kwargs)

    genomic_min = max(gene_to_draw.start - buffer_length, 1)
    genomic_max = gene_to_draw.end + buffer_length

    plots = []
    for axis_name, bam_file, density, ymax, stranded in read_depth_plots:
        # one plot per bam
        plots.append(
            bam_to_scatter(
                bam_file,
                gene_to_draw.chr,
                genomic_min,
                genomic_max,
                strand=gene_to_draw.get_strand() if stranded else None,
                ymax=ymax,
                density=density,
                axis_name=axis_name,
                min_mapping_quality=min_mapping_quality,
                ymax_color=ymax_color,
            )
        )

    for i, (marker_name, marker_start, marker_end) in enumerate(markers):
        markers[i] = _annotate.base.BioInterval(
            gene_to_draw.chr, marker_start, marker_end, name=marker_name
        )

    canvas = None
    attempts = 1
    while True:
        try:
            canvas = draw_multi_transcript_overlay(
                settings,
                gene_to_draw,
                vmarkers=markers,
                plots=plots,
                window_buffer=buffer_length,
                log=_util.LOG,
            )
            break
        except DrawingFitError as err:
            if attempts > max_drawing_retries:
                raise err
            _util.LOG('Drawing fit: extending window', drawing_width_iter_increase)
            settings.width += drawing_width_iter_increase
            attempts += 1

    svg_output_file = os.path.join(output, '{}_{}_overlay.svg'.format(gene_to_draw.name, gene_name))
    _util.LOG('writing:', svg_output_file)

    canvas.saveas(svg_output_file)


def convert_main(inputs, outputfile, file_type, strand_specific=False, assume_no_untemplated=True):
    bpp_results = convert_tool_output(
        inputs,
        file_type,
        strand_specific,
        _util.LOG,
        True,
        assume_no_untemplated=assume_no_untemplated,
    )
    if os.path.dirname(outputfile):
        _util.mkdirp(os.path.dirname(outputfile))
    _util.output_tabbed_file(bpp_results, outputfile)


def create_parser(argv):
    parser = argparse.ArgumentParser(formatter_class=_config.CustomHelpFormatter)
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version='%(prog)s version ' + __version__,
        help='Outputs the version number',
    )
    subp = parser.add_subparsers(
        dest='command', help='specifies which step/stage in the pipeline or which subprogram to use'
    )
    subp.required = True
    required = {}  # hold required argument group by subparser command name
    optional = {}  # hold optional argument group by subparser command name
    for command in SUBCOMMAND.values():
        subparser = subp.add_parser(
            command, formatter_class=_config.CustomHelpFormatter, add_help=False
        )
        required[command] = subparser.add_argument_group('required arguments')
        optional[command] = subparser.add_argument_group('optional arguments')
        optional[command].add_argument(
            '-h', '--help', action='help', help='show this help message and exit'
        )
        optional[command].add_argument(
            '-v',
            '--version',
            action='version',
            version='%(prog)s version ' + __version__,
            help='Outputs the version number',
        )
        optional[command].add_argument('--log', help='redirect stdout to a log file', default=None)
        optional[command].add_argument(
            '--log_level',
            help='level of logging to output',
            choices=['INFO', 'DEBUG'],
            default='INFO',
        )
        if command != SUBCOMMAND.CONVERT:
            optional[command].add_argument(
                '--config', '-c', help='path to the JSON config file', type=filepath, required=True
            )

    # convert
    required[SUBCOMMAND.CONVERT].add_argument(
        '--file_type',
        choices=sorted(SUPPORTED_TOOL.values()),
        required=True,
        help='Indicates the input file type to be parsed',
    )
    optional[SUBCOMMAND.CONVERT].add_argument(
        '--strand_specific', type=tab.cast_boolean, default=False
    )
    optional[SUBCOMMAND.CONVERT].add_argument(
        '--assume_no_untemplated', type=tab.cast_boolean, default=True
    )
    required[SUBCOMMAND.CONVERT].add_argument(
        '--outputfile', '-o', required=True, help='path to the outputfile', metavar='FILEPATH'
    )

    for command in set(SUBCOMMAND.values()) - {SUBCOMMAND.CONVERT}:
        required[command].add_argument(
            '-o', '--output', help='path to the output directory', required=True
        )

    # add the inputs argument
    for command in [
        SUBCOMMAND.CLUSTER,
        SUBCOMMAND.ANNOTATE,
        SUBCOMMAND.VALIDATE,
        SUBCOMMAND.PAIR,
        SUBCOMMAND.SUMMARY,
        SUBCOMMAND.CONVERT,
    ]:
        required[command].add_argument(
            '-n',
            '--inputs',
            nargs='+',
            help='path to the input files',
            required=True,
            metavar='FILEPATH',
        )

    # library specific commands
    for command in [SUBCOMMAND.CLUSTER, SUBCOMMAND.VALIDATE, SUBCOMMAND.ANNOTATE]:
        required[command].add_argument(
            '--library', '-l', required=True, help='The library to run the current step on'
        )

    # overlay arguments
    required[SUBCOMMAND.OVERLAY].add_argument('gene_name', help='Gene ID or gene alias to be drawn')
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--buffer_length',
        default=0,
        type=int,
        help='minimum genomic length to plot on either side of the target gene',
    )
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--read_depth_plot',
        dest='read_depth_plots',
        metavar='<axis name STR> <bam FILEPATH> [density FLOAT] [ymax INT] [stranded BOOL]',
        nmin=2,
        nmax=5,
        help='bam file to use as data for plotting read_depth',
        action=_config.RangeAppendAction,
    )
    optional[SUBCOMMAND.OVERLAY].add_argument(
        '--marker',
        dest='markers',
        metavar='<label STR> <start INT> [end INT]',
        nmin=2,
        nmax=3,
        help='Marker on the diagram given by genomic position, May be a single position or a range. '
        'The label should be a short descriptor to avoid overlapping labels on the diagram',
        action=_config.RangeAppendAction,
    )

    return parser, _util.MavisNamespace(**parser.parse_args(argv).__dict__)


def main(argv=None):
    """
    sets up the parser and checks the validity of command line args
    loads reference files and redirects into subcommand main functions

    Args:
        argv (list): List of arguments, defaults to command line arguments
    """
    if argv is None:  # need to do at run time or patching will not behave as expected
        argv = sys.argv[1:]
    start_time = int(time.time())
    parser, args = create_parser(argv)

    if args.command == SUBCOMMAND.OVERLAY:
        args = check_overlay_args(args, parser)

    log_conf = {'format': '{message}', 'style': '{', 'level': args.log_level}

    original_logging_handlers = logging.root.handlers[:]
    for handler in logging.root.handlers:
        logging.root.removeHandler(handler)
    if args.log:  # redirect stdout AND stderr to a log file
        log_conf['filename'] = args.log
    logging.basicConfig(**log_conf)

    _util.LOG('MAVIS: {}'.format(__version__))
    _util.LOG('hostname:', platform.node(), time_stamp=False)
    _util.log_arguments(args)

    config: Dict = dict()
    try:
        with open(args.config, 'r') as fh:
            config = json.load(fh)
            _config.validate_config(config, args.command)
    except AttributeError as err:
        if args.command != SUBCOMMAND.CONVERT:
            raise err

    if args.command == SUBCOMMAND.VALIDATE:
        args.aligner_version = get_aligner_version(config['validate']['aligner'])
    # try checking the input files exist
    try:
        args.inputs = _util.bash_expands(*args.inputs)
    except AttributeError:
        pass
    except FileNotFoundError:
        parser.error('--inputs file(s) for {} {} do not exist'.format(args.command, args.inputs))

    # decide which main function to execute
    command = args.command

    try:
        if command == SUBCOMMAND.CLUSTER:
            cluster_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
                library=args.library,
            )
        elif command == SUBCOMMAND.VALIDATE:
            validate_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
                library=args.library,
            )
        elif command == SUBCOMMAND.ANNOTATE:
            annotate_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
                library=args.library,
            )
        elif command == SUBCOMMAND.PAIR:
            pairing_main.main(
                inputs=args.inputs,
                output=args.output,
                start_time=start_time,
                config=config,
            )
        elif command == SUBCOMMAND.SUMMARY:
            summary_main.main(**args, start_time=start_time)
        elif command == SUBCOMMAND.CONVERT:
            convert_main(
                args.inputs,
                args.outputfile,
                args.file_type,
                args.strand_specific,
                args.assume_no_untemplated,
            )
        else:
            overlay_main(**args)

        duration = int(time.time()) - start_time
        hours = duration - duration % 3600
        minutes = duration - hours - (duration - hours) % 60
        seconds = duration - hours - minutes
        _util.LOG(
            'run time (hh/mm/ss): {}:{:02d}:{:02d}'.format(hours // 3600, minutes // 60, seconds),
            time_stamp=False,
        )
        _util.LOG('run time (s): {}'.format(duration), time_stamp=False)
    except Exception as err:
        raise err
    finally:
        try:
            for handler in logging.root.handlers:
                logging.root.removeHandler(handler)
            for handler in original_logging_handlers:
                logging.root.addHandler(handler)
        except Exception as err:
            print(err)


if __name__ == '__main__':
    main()
