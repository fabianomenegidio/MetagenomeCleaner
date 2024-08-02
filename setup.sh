#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Create the directory structure
mkdir -p metagenome_cleaner
mkdir -p metagenome_cleaner/scripts
mkdir -p metagenome_cleaner/databases/bowtie2
mkdir -p metagenome_cleaner/databases/kraken2
mkdir -p metagenome_cleaner/logs
mkdir -p metagenome_cleaner/reports
mkdir -p metagenome_cleaner/output
mkdir -p metagenome_cleaner/tests
mkdir -p metagenome_cleaner/docs

# Create a README file
cat << 'EOF' > metagenome_cleaner/README.md
# Metagenome Cleaner

Metagenome Cleaner is a tool designed to clean metagenome FASTQ libraries by integrating various modules.

## Features
- Quality filtering
- Contaminant removal (PhiX, Human, Host)
- Parallel execution of modules
- Comprehensive reports with MultiQC

## Installation
- Ensure that the required tools are installed: `fastp`, `kraken2`, `bowtie2`, `parallel`, `multiqc`, `fastqc`
- To install these tools, you can use conda or your package manager.

## Usage
Run the tool from the command line and follow the menu options.

## Development
This tool was developed by the Laboratory of Bioinformatics and Omic Sciences (LaBiOmicS) at Universidade de Mogi das Cruzes (UMC), in partnership with the Universidade Federal de Viçosa (UFV), Universidade Federal do ABC (UFABC), and Universidade Estadual do Sudoeste da Bahia (UESB).
EOF

# Create a setup script for Python packages
cat << 'EOF' > metagenome_cleaner/setup.py
from setuptools import setup, find_packages

setup(
    name='metagenome_cleaner',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'argparse',
        'configparser',
        'datetime',
        'glob2',
        'shutil',
        'subprocess'
    ],
    entry_points={
        'console_scripts': [
            'metagenome_cleaner=metagenome_cleaner:main',
        ],
    },
)
EOF

# Create a requirements file
cat << 'EOF' > metagenome_cleaner/requirements.txt
argparse
configparser
datetime
glob2
shutil
subprocess
EOF

# Create a configuration file
cat << 'EOF' > metagenome_cleaner/config.ini
[databases]
bowtie2_phix_index=metagenome_cleaner/databases/bowtie2/phix_index
bowtie2_human_index=metagenome_cleaner/databases/bowtie2/human_index
kraken2_human_db=metagenome_cleaner/databases/kraken2/human_db
bowtie2_dog_index=metagenome_cleaner/databases/bowtie2/dog_index
kraken2_dog_db=metagenome_cleaner/databases/kraken2/dog_db
bowtie2_cat_index=metagenome_cleaner/databases/bowtie2/cat_index
kraken2_cat_db=metagenome_cleaner/databases/kraken2/cat_db
bowtie2_rat_index=metagenome_cleaner/databases/bowtie2/rat_index
kraken2_rat_db=metagenome_cleaner/databases/kraken2/rat_db
bowtie2_mouse_index=metagenome_cleaner/databases/bowtie2/mouse_index
kraken2_mouse_db=metagenome_cleaner/databases/kraken2/mouse_db
bowtie2_cow_index=metagenome_cleaner/databases/bowtie2/cow_index
kraken2_cow_db=metagenome_cleaner/databases/kraken2/cow_db
bowtie2_pig_index=metagenome_cleaner/databases/bowtie2/pig_index
kraken2_pig_db=metagenome_cleaner/databases/kraken2/pig_db
bowtie2_horse_index=metagenome_cleaner/databases/bowtie2/horse_index
kraken2_horse_db=metagenome_cleaner/databases/kraken2/horse_db
bowtie2_zebrafish_index=metagenome_cleaner/databases/bowtie2/zebrafish_index
kraken2_zebrafish_db=metagenome_cleaner/databases/kraken2/zebrafish_db
bowtie2_yeast_index=metagenome_cleaner/databases/bowtie2/yeast_index
kraken2_yeast_db=metagenome_cleaner/databases/kraken2/yeast_db
EOF

# Create the main Python script
cat << 'EOF' > metagenome_cleaner/metagenome_cleaner.py
import argparse
import subprocess
import shutil
import os
import glob
import configparser
from datetime import datetime
from collections import defaultdict

def read_config(config_file='config.ini'):
    """
    Read the configuration file.
    """
    config = configparser.ConfigParser()
    config.read(config_file)
    return config

def main():
    """
    Main function to handle the command-line interface and execute the appropriate functions.
    """
    print("Welcome to Metagenome Cleaner!")
    print("Please select an option from the menu below:")
    print("1. About")
    print("2. Quality Filtering")
    print("3. Contaminant Removal")
    print("4. Parallel Execution (Optional)")
    print("5. Help")
    print("\nUse -h or --help for more information on each option.\n")

    parser = argparse.ArgumentParser(description='Metagenome Cleaner - A tool for cleaning metagenome FASTQ libraries.')
    subparsers = parser.add_subparsers(dest='command')

    # About command
    about_parser = subparsers.add_parser('about', help='Display information about Metagenome Cleaner')
    about_parser.set_defaults(func=display_about)

    # Quality Filtering command
    quality_filtering_parser = subparsers.add_parser('quality', help='Perform quality filtering on FASTQ files')
    quality_filtering_parser.add_argument('-1', '--forward', type=str, help='Specify the forward FASTQ file for paired-end or single-end input')
    quality_filtering_parser.add.argument('-2', '--reverse', type=str, help='Specify the reverse FASTQ file for paired-end input')
    quality_filtering_parser.add.argument('-U', '--single', type=str, help='Specify the single-end FASTQ file')
    quality_filtering_parser.add.argument('--input-dir', type=str, help='Specify the input directory containing FASTQ files')
    quality_filtering_parser.add.argument('--output', type=str, required=True, help='Specify the output directory')
    quality_filtering_parser.add.argument('--save-intermediates', action='store_true', help='Save all intermediate files')

    # Add all relevant fastp parameters
    quality_filtering_parser.add.argument('--qualified_quality_phred', type=int, default=20, help='Quality value that a base is qualified (default: 20) (fastp)')
    quality_filtering_parser.add.argument('--unqualified_percent_limit', type=float, help='How many percents of unqualified bases are allowed in a read (fastp)')
    quality_filtering_parser.add.argument('--n_base_limit', type=int, help='Maximum number of N bases allowed in a read (fastp)')
    quality_filtering_parser.add.argument('--length_required', type=int, help='Reads shorter than this length will be discarded (fastp)')
    quality_filtering_parser.add.argument('--adapter_sequence', type=str, help='Adapter sequence for read1 (fastp)')
    quality_filtering_parser.add.argument('--adapter_sequence_r2', type=str, help='Adapter sequence for read2 (paired-end) (fastp)')
    quality_filtering_parser.add.argument('--trim_poly_g', action='store_true', help='Enable polyG tail trimming (useful for NovaSeq data) (fastp)')
    quality_filtering_parser.add.argument('--trim_poly_x', action='store_true', help='Enable polyX tail trimming (fastp)')
    quality_filtering_parser.add.argument('--cut_front', action='store_true', help='Enable front end trimming (fastp)')
    quality_filtering_parser.add.argument('--cut_tail', action='store_true', help='Enable tail end trimming (fastp)')
    quality_filtering_parser.add.argument('--cut_window_size', type=int, help='Window size for sliding window cutting (fastp)')
    quality_filtering_parser.add.argument('--cut_mean_quality', type=int, help='Mean quality requirement for sliding window cutting (fastp)')
    quality_filtering_parser.add.argument('--html', type=str, help='File name to store HTML report (fastp)')
    quality_filtering_parser.add.argument('--json', type=str, help='File name to store JSON report (fastp)')
    quality_filtering_parser.add.argument('--thread', type=int, help='Number of worker threads to be used (fastp)')

    quality_filtering_parser.set_defaults(func=quality_filtering)

    # Contaminant Removal command
    contaminant_removal_parser = subparsers.add_parser('contaminant_removal', help='Remove contaminant organisms from FASTQ files')
    contaminant_removal_subparsers = contaminant_removal_parser.add_subparsers(dest='contaminant_type')

    # Submenu to remove PhiX
    phix_parser = contaminant_removal_subparsers.add_parser('phix', help='Remove PhiX contaminant sequences')
    phix_parser.add.argument('-1', '--forward', type=str, help='Specify the forward FASTQ file for paired-end or single-end input')
    phix_parser.add.argument('-2', '--reverse', type=str, help='Specify the reverse FASTQ file for paired-end input')
    phix_parser.add.argument('-U', '--single', type=str, help='Specify the single-end FASTQ file')
    phix_parser.add.argument('--input-dir', type=str, help='Specify the input directory containing FASTQ files')
    phix_parser.add.argument('--output', type=str, required=True, help='Specify the output directory')
    phix_parser.add.argument('--save-intermediates', action='store_true', help='Save all intermediate files')
    phix_parser.add.argument('--N', type=int, help='Number of mismatches allowed in a seed alignment (bowtie2)')
    phix_parser.add.argument('--L', type=int, help='Length of seed substrings (bowtie2)')
    phix_parser.add.argument('--i', type=str, help='Interval function for seed extension (bowtie2)')
    phix_parser.add.argument('--rdg', type=str, help='Read gap open and extend penalties (bowtie2)')
    phix_parser.add.argument('--rfg', type=str, help='Reference gap open and extend penalties (bowtie2)')
    phix_parser.add.argument('--score-min', type=str, help='Minimum alignment score (bowtie2)')
    phix_parser.add.argument('--very-fast', action='store_true', help='Use very fast preset (bowtie2)')
    phix_parser.add.argument('--fast', action='store_true', help='Use fast preset (bowtie2)')
    phix_parser.add.argument('--sensitive', action='store_true', help='Use sensitive preset (bowtie2)')
    phix_parser.add.argument('--very-sensitive', action='store_true', help='Use very sensitive preset (bowtie2)')
    phix_parser.add.argument('--al', type=str, help='File to write aligned reads (single-end) (bowtie2)')
    phix_parser.add.argument('--al-conc', type=str, help='Files to write aligned reads (paired-end) (bowtie2)')
    phix_parser.add.argument('--un', type=str, help='File to write unaligned reads (single-end) (bowtie2)')
    phix_parser.add.argument('--un-conc', type=str, help='Files to write unaligned reads (paired-end) (bowtie2)')
    phix_parser.add.argument('--k', type=int, help='Report up to <int> valid alignments per read (bowtie2)')
    phix_parser.add.argument('--a', action='store_true', help='Report all alignments (bowtie2)')
    phix_parser.add.argument('--maxins', type=int, help='Maximum fragment length for paired-end reads (bowtie2)')
    phix_parser.set_defaults(func=remove_phix_contaminants)

    # Submenu to remove human contaminants
    human_parser = contaminant_removal_subparsers.add_parser('human', help='Remove Human contaminant sequences')
    human_parser.add.argument('-1', '--forward', type=str, help='Specify the forward FASTQ file for paired-end or single-end input')
    human_parser.add.argument('-2', '--reverse', type=str, help='Specify the reverse FASTQ file for paired-end input')
    human_parser.add.argument('-U', '--single', type=str, help='Specify the single-end FASTQ file')
    human_parser.add.argument('--input-dir', type=str, help='Specify the input directory containing FASTQ files')
    human_parser.add.argument('--output', type=str, required=True, help='Specify the output directory')
    human_parser.add.argument('--save-intermediates', action='store_true', help='Save all intermediate files')

    # Classification and alignment parameters
    human_parser.add.argument('--threads', type=int, help='Number of threads (bowtie2, kraken2)')
    human_parser.add.argument('--N', type=int, help='Number of mismatches allowed in a seed alignment (bowtie2)')
    human_parser.add.argument('--L', type=int, help='Length of seed substrings (bowtie2)')
    human_parser.add.argument('--i', type=str, help='Interval function for seed extension (bowtie2)')
    human_parser.add.argument('--rdg', type=str, help='Read gap open and extend penalties (bowtie2)')
    human_parser.add.argument('--rfg', type=str, help='Reference gap open and extend penalties (bowtie2)')
    human_parser.add.argument('--score-min', type=str, help='Minimum alignment score (bowtie2)')
    human_parser.add.argument('--very-fast', action='store_true', help='Use very fast preset (bowtie2)')
    human_parser.add.argument('--fast', action='store_true', help='Use fast preset (bowtie2)')
    human_parser.add.argument('--sensitive', action='store_true', help='Use sensitive preset (bowtie2)')
    human_parser.add.argument('--very-sensitive', action='store_true', help='Use very sensitive preset (bowtie2)')
    human_parser.add.argument('--al', type=str, help='File to write aligned reads (single-end) (bowtie2)')
    human_parser.add.argument('--al-conc', type=str, help='Files to write aligned reads (paired-end) (bowtie2)')
    human_parser.add.argument('--un', type=str, help='File to write unaligned reads (single-end) (bowtie2)')
    human_parser.add.argument('--un-conc', type=str, help='Files to write unaligned reads (paired-end) (bowtie2)')
    human_parser.add.argument('--k', type=int, help='Report up to <int> valid alignments per read (bowtie2)')
    human_parser.add.argument('--a', action='store_true', help='Report all alignments (bowtie2)')
    human_parser.add.argument('--maxins', type=int, help='Maximum fragment length for paired-end reads (bowtie2)')

    # Kraken2 parameters
    human_parser.add.argument('--confidence', type=float, help='Confidence score threshold (kraken2)')
    human_parser.add.argument('--minimum-hit-groups', type=int, help='Minimum number of hit groups (kraken2)')
    human_parser.add.argument('--report', type=str, help='File to write report (kraken2)')
    human_parser.add.argument('--use-names', action='store_true', help='Use scientific names in report (kraken2)')
    human_parser.add.argument('--minimum-base-quality', type=int, help='Minimum base quality (kraken2)')
    human_parser.add.argument('--report-minimizer-data', action='store_true', help='Include minimizer data in report (kraken2)')

    human_parser.set_defaults(func=remove_human_contaminants)

    # Submenu to remove host contaminants
    host_parser = contaminant_removal_subparsers.add_parser('host', help='Remove host contaminant sequences')
    host_parser.add.argument('--host', type=str, required=True, choices=['dog', 'cat', 'rat', 'mouse', 'cow', 'pig', 'horse', 'zebrafish', 'yeast', 'custom'], help='Specify the host organism (dog, cat, rat, mouse, cow, pig, horse, zebrafish, yeast, custom)')
    host_parser.add.argument('--db-path-bowtie2', type=str, help='Path to the custom Bowtie2 database (Required if --host is custom)')
    host_parser.add.argument('--db-path-kraken2', type=str, help='Path to the custom Kraken2 database (Required if --host is custom)')
    host_parser.add.argument('-1', '--forward', type=str, help='Specify the forward FASTQ file for paired-end or single-end input')
    host_parser.add.argument('-2', '--reverse', type=str, help='Specify the reverse FASTQ file for paired-end input')
    host_parser.add.argument('-U', '--single', type.str, help='Specify the single-end FASTQ file')
    host_parser.add.argument('--input-dir', type.str, help='Specify the input directory containing FASTQ files')
    host_parser.add.argument('--output', type.str, required=True, help='Specify the output directory')
    host_parser.add.argument('--save-intermediates', action='store_true', help='Save all intermediate files')

    # Classification and alignment parameters
    host_parser.add.argument('--threads', type.int, help='Number of threads (bowtie2, kraken2)')
    host_parser.add.argument('--N', type.int, help='Number of mismatches allowed in a seed alignment (bowtie2)')
    host_parser.add.argument('--L', type.int, help='Length of seed substrings (bowtie2)')
    host_parser.add.argument('--i', type.str, help='Interval function for seed extension (bowtie2)')
    host_parser.add.argument('--rdg', type.str, help='Read gap open and extend penalties (bowtie2)')
    host_parser.add.argument('--rfg', type.str, help='Reference gap open and extend penalties (bowtie2)')
    host_parser.add.argument('--score-min', type.str, help='Minimum alignment score (bowtie2)')
    host_parser.add.argument('--very-fast', action='store_true', help='Use very fast preset (bowtie2)')
    host_parser.add.argument('--fast', action='store_true', help='Use fast preset (bowtie2)')
    host_parser.add.argument('--sensitive', action='store_true', help='Use sensitive preset (bowtie2)')
    host_parser.add.argument('--very-sensitive', action='store_true', help='Use very sensitive preset (bowtie2)')
    host_parser.add.argument('--al', type.str, help='File to write aligned reads (single-end) (bowtie2)')
    host_parser.add.argument('--al-conc', type.str, help='Files to write aligned reads (paired-end) (bowtie2)')
    host_parser.add.argument('--un', type.str, help='File to write unaligned reads (single-end) (bowtie2)')
    host_parser.add.argument('--un-conc', type.str, help='Files to write unaligned reads (paired-end) (bowtie2)')
    host_parser.add.argument('--k', type.int, help='Report up to <int> valid alignments per read (bowtie2)')
    host_parser.add.argument('--a', action='store_true', help='Report all alignments (bowtie2)')
    host_parser.add.argument('--maxins', type.int, help='Maximum fragment length for paired-end reads (bowtie2)')

    # Kraken2 parameters
    host_parser.add.argument('--confidence', type.float, help='Confidence score threshold (kraken2)')
    host_parser.add.argument('--minimum-hit-groups', type.int, help='Minimum number of hit groups (kraken2)')
    host_parser.add.argument('--report', type.str, help='File to write report (kraken2)')
    host_parser.add.argument('--use-names', action='store_true', help='Use scientific names in report (kraken2)')
    host_parser.add.argument('--minimum-base-quality', type.int, help='Minimum base quality (kraken2)')
    host_parser.add.argument('--report-minimizer-data', action='store_true', help='Include minimizer data in report (kraken2)')

    host_parser.set_defaults(func=remove_host_contaminants)

    # Parallel command
    parallel_parser = subparsers.add_parser('parallel', help='Execute modules in parallel (Optional)')
    parallel_subparsers = parallel_parser.add_subparsers(dest='parallel_command')

    # Sub-menu for parallel execution
    parallel_exec_parser = parallel_subparsers.add_parser('exec', help='Execute a module in parallel')
    parallel_exec_parser.add.argument('--module', type=str, required=True, choices=['quality', 'contaminant_removal'], help='Specify the module to execute in parallel (Available: quality, contaminant_removal)')
    parallel_exec_parser.add.argument('--contaminant_type', type.str, choices=['phix', 'human', 'host'], help='Specify the type of contaminant to remove (phix, human, host)')
    parallel_exec_parser.add.argument('--host', type.str, choices=['dog', 'cat', 'rat', 'mouse', 'cow', 'pig', 'horse', 'zebrafish', 'yeast', 'custom'], help='Specify the host organism for host contaminant removal (dog, cat, rat, mouse, cow, pig, horse, zebrafish, yeast, custom)')
    parallel_exec_parser.add.argument('--db-path-bowtie2', type.str, help='Path to the custom Bowtie2 database (Required if --host is custom)')
    parallel_exec_parser.add.argument('--db-path-kraken2', type.str, help='Path to the custom Kraken2 database (Required if --host is custom)')
    parallel_exec_parser.add.argument('--params', type.str, help='Specify any additional parameters for parallel execution')
    parallel_exec_parser.set_defaults(func=execute_parallel)

    args = parser.parse_args()
    if args.command:
        args.func(args)
    else:
        parser.print_help()

def display_about(args):
    """
    Display information about the tool.
    """
    print("Metagenome Cleaner is a tool designed to clean metagenome FASTQ libraries by integrating various modules.")
    print("It offers parallel execution of modules for efficient data processing.")
    print("\nThis tool was developed by the Laboratory of Bioinformatics and Omic Sciences (LaBiOmicS) at Universidade de Mogi das Cruzes (UMC),")
    print("in partnership with the Universidade Federal de Viçosa (UFV), Universidade Federal do ABC (UFABC), and Universidade Estadual do Sudoeste da Bahia (UESB).")

def display_help(args):
    """
    Display help information about the tool.
    """
    print("Metagenome Cleaner Help:")
    print("This tool helps in cleaning metagenome FASTQ libraries.")
    print("Use the 'quality' option to perform quality filtering on FASTQ files.")
    print("Use the 'contaminant_removal' option to remove contaminants from FASTQ files.")
    print("Use the 'parallel' option to execute modules in parallel (Optional).")
    print("For more information on each option, use the -h or --help flag.")

def check_tool_installed(tool_name):
    """
    Check if a tool is installed on the system.
    """
    return shutil.which(tool_name) is not None

def suggest_installation(tool_name):
    """
    Suggest installation command for a missing tool.
    """
    print(f"{tool_name} is not installed. Please install it to proceed.")
    print(f"To install {tool_name}, you can use the following command:")
    if tool_name == 'fastp':
        print("  conda install -c bioconda fastp")
    elif tool_name == 'kraken2':
        print("  conda install -c bioconda kraken2")
    elif tool_name == 'bowtie2':
        print("  conda install -c bioconda bowtie2")
    elif tool_name == 'parallel':
        print("  sudo apt-get install parallel")
    elif tool_name == 'multiqc':
        print("  conda install -c bioconda multiqc")
    elif tool_name == 'fastqc':
        print("  conda install -c bioconda fastqc")

def create_project_dir():
    """
    Create a project directory with a timestamp.
    """
    timestamp = datetime.now().strftime("%d%m%y%H%M")
    project_dir = f"Project_{timestamp}"
    os.makedirs(project_dir, exist_ok=True)
    return project_dir

def create_module_dir(project_dir, module_name):
    """
    Create a module directory within the project directory.
    """
    module_dir = os.path.join(project_dir, module_name)
    os.makedirs(module_dir, exist_ok=True)
    return module_dir

def run_fastqc(input_files, output_dir):
    """
    Run FastQC on the given input files and save the reports to the output directory.
    """
    if not check_tool_installed('fastqc'):
        suggest_installation('fastqc')
        return

    commands = [f"fastqc -o {output_dir} {input_file}" for input_file in input_files]
    run_commands(commands)

def quality_filtering(args):
    """
    Perform quality filtering on FASTQ files using fastp.
    """
    if not check_tool_installed('fastp'):
        suggest_installation('fastp')
        return

    forward = args.forward
    reverse = args.reverse
    single = args.single
    input_dir = args.input_dir
    output_dir = args.output
    save_intermediates = args.save_intermediates

    params = []

    project_dir = create_project_dir()
    module_dir = create_module_dir(project_dir, 'quality')

    # Run FastQC before filtering
    input_files = []
    if input_dir:
        input_files = glob.glob(os.path.join(input_dir, '*.fastq*'))
    else:
        if forward:
            input_files.append(forward)
        if reverse:
            input_files.append(reverse)
        if single:
            input_files.append(single)
    run_fastqc(input_files, module_dir)

    if args.qualified_quality_phred:
        params.append(f("--qualified_quality_phred {args.qualified_quality_phred}")
    if args.unqualified_percent_limit:
        params.append(f("--unqualified_percent_limit {args.unqualified_percent_limit}")
    if args.n_base_limit:
        params.append(f("--n_base_limit {args.n_base_limit}")
    if args.length_required:
        params.append(f("--length_required {args.length_required}")
    if args.adapter_sequence:
        params.append(f("--adapter_sequence {args.adapter_sequence}")
    if args.adapter_sequence_r2:
        params.append(f("--adapter_sequence_r2 {args.adapter_sequence_r2}")
    if args.trim_poly_g:
        params.append("--trim_poly_g")
    if args.trim_poly_x:
        params.append("--trim_poly_x")
    if args.cut_front:
        params.append("--cut_front")
    if args.cut_tail:
        params.append("--cut_tail")
    if args.cut_window_size:
        params.append(f("--cut_window_size {args.cut_window_size}")
    if args.cut_mean_quality:
        params.append(f("--cut_mean_quality {args.cut_mean_quality}")
    if args.html:
        params.append(f("--html {os.path.join(module_dir, args.html)}")
    if args.json:
        params.append(f("--json {os.path.join(module_dir, args.json)}")
    if args.thread:
        params.append(f("--thread {args.thread}")

    params_str = " ".join(params)

    commands = []

    def get_output_filename(input_file, suffix):
        base_name = os.path.basename(input_file).replace('.fastq', '').replace('.gz', '').replace('_R1', '').replace('_R2', '').replace('_1', '').replace('_2', '')
        return os.path.join(module_dir, f"{base_name}_filtered_{suffix}.fastq")

    if input_dir:
        paired_files = defaultdict(list)
        single_files = []

        files = glob.glob(os.path.join(input_dir, '*.fastq*'))
        for file in files:
            if '_R1' in file or '_1' in file:
                paired_files[file.replace('_R1', '').replace('_1', '')).append(file)
            elif '_R2' in file or '_2' in file:
                paired_files[file.replace('_R2', '').replace('_2', '')).append(file)
            else:
                single_files.append(file)

        for base, pair in paired_files.items():
            if len(pair) == 2:
                out1 = get_output_filename(pair[0], "1")
                out2 = get_output_filename(pair[1], "2")
                command = f"fastp -i {pair[0]} -I {pair[1]} -o {out1} -O {out2} {params_str}"
                commands.append(command)

        for file in single_files:
            out = get_output_filename(file, "")
            command = f"fastp -i {file} -o {out} {params_str}"
            commands.append(command)

    else:
        if forward and reverse:
            out1 = get_output_filename(forward, "1")
            out2 = get_output_filename(reverse, "2")
            command = f"fastp -i {forward} -I {reverse} -o {out1} -O {out2} {params_str}"
        elif single:
            out = get_output_filename(single, "")
            command = f"fastp -i {single} -o {out} {params_str}"
        else:
            print("Error: Please provide either --forward and --reverse for paired-end data, --single for single-end data, or --input-dir for a directory of FASTQ files.")
            return
        commands.append(command)

    run_commands(commands, module_dir, save_intermediates)

    # Run FastQC after filtering
    filtered_files = glob.glob(os.path.join(module_dir, '*_filtered*.fastq'))
    run_fastqc(filtered_files, module_dir)

    # Generate MultiQC report
    generate_multiqc_report(module_dir)

def run_commands(commands, module_dir=None, save_intermediates=False):
    """
    Execute a list of shell commands and log their execution.
    """
    log_file = os.path.join(module_dir, "execution_log.txt") if module_dir else "execution_log.txt"
    with open(log_file, 'a') as log:
        for command in commands:
            print(f"Executing command: {command}")
            log.write(f"Executing command: {command}\n")
            try:
                subprocess.run(command, shell=True, check=True)
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while executing the command: {e}")
                log.write(f"An error occurred while executing the command: {e}\n")

    if module_dir and not save_intermediates:
        # Remove intermediate files if not saving intermediates
        for file in os.listdir(module_dir):
            if not file.endswith('_filtered.fastq') and not file.endswith('_log.txt') and not file.endswith('.html') and not file.endswith('.json') and not file.endswith('_fastqc.zip'):
                os.remove(os.path.join(module_dir, file))

def execute_parallel(args):
    """
    Execute modules in parallel using GNU Parallel.
    """
    if not check_tool_installed('parallel'):
        suggest_installation('parallel')
        return

    module = args.module
    params = args.params if args.params else ""

    if module == "quality":
        command = f"parallel {params} ::: 'metagenome_cleaner quality --forward {{}} --reverse {{}} --output {{}}_out'"
    elif module == "contaminant_removal":
        contaminant_type = args.contaminant_type
        if contaminant_type == 'host':
            host = args.host
            if host == 'custom':
                db_path_bowtie2 = args.db_path_bowtie2
                db_path_kraken2 = args.db_path_kraken2
                command = f"parallel {params} ::: 'metagenome_cleaner contaminant_removal host --host custom --db-path-bowtie2 {db_path_bowtie2} --db-path-kraken2 {db_path_kraken2} --forward {{}} --reverse {{}} --output {{}}_out'"
            else:
                command = f"parallel {params} ::: 'metagenome_cleaner contaminant_removal host --host {host} --forward {{}} --reverse {{}} --output {{}}_out'"
        else:
            command = f"parallel {params} ::: 'metagenome_cleaner contaminant_removal {contaminant_type} --forward {{}} --reverse {{}} --output {{}}_out'"
    else:
        print("Error: Unsupported module for parallel execution.")
        return

    print(f"Executing module in parallel with command: {command}")
    try:
        subprocess.run(command, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while executing the command: {e}")

def remove_phix_contaminants(args):
    """
    Remove PhiX contaminant sequences using Bowtie2.
    """
    config = read_config()
    phix_index = config['databases']['bowtie2_phix_index']

    if not check_tool_installed('bowtie2'):
        suggest_installation('bowtie2')
        return

    forward = args.forward
    reverse = args.reverse
    single = args.single
    input_dir = args.input_dir
    output_dir = args.output
    save_intermediates = args.save_intermediates
    align_params = ["--very-sensitive-local"]  # Default parameter

    if args.N:
        align_params.append(f("-N {args.N}")
    if args.L:
        align_params.append(f("-L {args.L}")
    if args.i:
        align_params.append(f("-i {args.i}")
    if args.rdg:
        align_params.append(f("--rdg {args.rdg}")
    if args.rfg:
        align_params.append(f("--rfg {args.rfg}")
    if args.score_min:
        align_params.append(f("--score-min {args.score_min}")
    if args.very_fast:
        align_params.append("--very-fast")
    if args.fast:
        align_params.append("--fast")
    if args.sensitive:
        align_params.append("--sensitive")
    if args.very_sensitive:
        align_params.append("--very-sensitive")
    if args.al:
        align_params.append(f("--al {args.al}")
    if args.al_conc:
        align_params.append(f("--al-conc {args.al_conc}")
    if args.un:
        align_params.append(f("--un {args.un}")
    if args.un_conc:
        align_params.append(f("--un-conc {args.un_conc}")
    if args.k:
        align_params.append(f("-k {args.k}")
    if args.a:
        align_params.append("-a")
    if args.maxins:
        align_params.append(f("--maxins {args.maxins}")

    align_params_str = " ".join(align_params)

    project_dir = create_project_dir()
    module_dir = create_module_dir(project_dir, 'phix_removal')

    commands = []

    def get_output_filename(input_file, suffix):
        base_name = os.path.basename(input_file).replace('.fastq', '').replace('.gz', '').replace('_R1', '').replace('_R2', '').replace('_1', '').replace('_2', '')
        return os.path.join(module_dir, f"{base_name}_phix_cleaned_{suffix}.fastq")

    if input_dir:
        paired_files = defaultdict(list)
        single_files = []

        files = glob.glob(os.path.join(input_dir, '*.fastq*'))
        for file in files:
            if '_R1' in file or '_1' in file:
                paired_files[file.replace('_R1', '').replace('_1', '')).append(file)
            elif '_R2' in file or '_2' in file:
                paired_files[file.replace('_R2', '').replace('_2', '')).append(file)
            else:
                single_files.append(file)

        for base, pair in paired_files.items():
            if len(pair) == 2:
                out1 = get_output_filename(pair[0], "1")
                out2 = get_output_filename(pair[1], "2")
                command = f"bowtie2 -x {phix_index} -1 {pair[0]} -2 {pair[1]} --un-conc {out1},{out2} {align_params_str}"
                commands.append(command)

        for file in single_files:
            out = get_output_filename(file, "")
            command = f"bowtie2 -x {phix_index} -U {file} --un {out} {align_params_str}"
            commands.append(command)

    else:
        if forward and reverse:
            out1 = get_output_filename(forward, "1")
            out2 = get_output_filename(reverse, "2")
            command = f"bowtie2 -x {phix_index} -1 {forward} -2 {reverse} --un-conc {out1},{out2} {align_params_str}"
        elif single:
            out = get_output_filename(single, "")
            command = f"bowtie2 -x {phix_index} -U {single} --un {out} {align_params_str}"
        else:
            print("Error: Please provide either --forward and --reverse for paired-end data, --single for single-end data, or --input-dir for a directory of FASTQ files.")
            return
        commands.append(command)

    run_commands(commands, module_dir, save_intermediates)

    # Run FastQC after removal
    cleaned_files = glob.glob(os.path.join(module_dir, '*_phix_cleaned*.fastq'))
    run_fastqc(cleaned_files, module_dir)

    # Generate MultiQC report
    generate_multiqc_report(module_dir)

def remove_human_contaminants(args):
    """
    Remove Human contaminant sequences using Kraken2 and Bowtie2.
    """
    config = read_config()
    human_db = config['databases']['kraken2_human_db']
    human_index = config['databases']['bowtie2_human_index']

    if not check_tool_installed('kraken2'):
        suggest_installation('kraken2')
        return
    if not check_tool_installed('bowtie2'):
        suggest_installation('bowtie2')
        return

    forward = args.forward
    reverse = args.reverse
    single = args.single
    input_dir = args.input_dir
    output_dir = args.output
    save_intermediates = args.save_intermediates
    threads = args.threads if args.threads else 1
    align_params = ["--very-sensitive-local"]  # Default parameter
    kraken2_params = []

    # Kraken2 parameters
    kraken2_params.append(f("--threads {threads}")
    if args.confidence:
        kraken2_params.append(f("--confidence {args.confidence}")
    if args.minimum_hit_groups:
        kraken2_params.append(f("--minimum-hit-groups {args.minimum_hit_groups}")
    if args.report:
        kraken2_params.append(f("--report {args.report}")
    if args.use_names:
        kraken2_params.append("--use-names")
    if args.minimum_base_quality:
        kraken2_params.append(f("--minimum-base-quality {args.minimum_base_quality}")
    if args.report_minimizer_data:
        kraken2_params.append("--report-minimizer-data")

    kraken2_params_str = " ".join(kraken2_params)

    # Bowtie2 parameters
    if args.N:
        align_params.append(f("-N {args.N}")
    if args.L:
        align_params.append(f("-L {args.L}")
    if args.i:
        align_params.append(f("-i {args.i}")
    if args.rdg:
        align_params.append(f("--rdg {args.rdg}")
    if args.rfg:
        align_params.append(f("--rfg {args.rfg}")
    if args.score_min:
        align_params.append(f("--score-min {args.score_min}")
    if args.very_fast:
        align_params.append("--very-fast")
    if args.fast:
        align_params.append("--fast")
    if args.sensitive:
        align_params.append("--sensitive")
    if args.very_sensitive:
        align_params.append("--very-sensitive")
    if args.al:
        align_params.append(f("--al {args.al}")
    if args.al_conc:
        align_params.append(f("--al-conc {args.al_conc}")
    if args.un:
        align_params.append(f("--un {args.un}")
    if args.un_conc:
        align_params.append(f("--un-conc {args.un_conc}")
    if args.k:
        align_params.append(f("-k {args.k}")
    if args.a:
        align_params.append("-a")
    if args.maxins:
        align_params.append(f("--maxins {args.maxins}")

    align_params_str = " ".join(align_params)

    project_dir = create_project_dir()
    module_dir = create_module_dir(project_dir, 'human_removal')

    commands = []

    def get_kraken2_output_filename(input_file):
        base_name = os.path.basename(input_file).replace('.fastq', '').replace('.gz', '')
        return os.path.join(module_dir, f"{base_name}_kraken2_unclassified.fastq")

    def get_bowtie2_output_filename(input_file, suffix):
        base_name = os.path.basename(input_file).replace('.fastq', '').replace('.gz', '').replace('_R1', '').replace('_R2', '').replace('_1', '').replace('_2', '')
        return os.path.join(module_dir, f"{base_name}_human_bowtie2_unaligned_{suffix}.fastq")

    def get_report_filename(input_file, suffix):
        base_name = os.path.basename(input_file).replace('.fastq', '').replace('.gz', '')
        return os.path.join(module_dir, f"{base_name}_{suffix}.txt")

    if input_dir:
        paired_files = defaultdict(list)
        single_files = []

        files = glob.glob(os.path.join(input_dir, '*.fastq*'))
        for file in files:
            if '_R1' in file or '_1' in file:
                paired_files[file.replace('_R1', '').replace('_1', '')).append(file)
            elif '_R2' in file or '_2' in file:
                paired_files[file.replace('_R2', '').replace('_2', '')).append(file)
            else:
                single_files.append(file)

        for base, pair in paired_files.items():
            if len(pair) == 2:
                kraken2_out1 = get_kraken2_output_filename(pair[0])
                kraken2_out2 = get_kraken2_output_filename(pair[1])
                kraken2_report = get_report_filename(pair[0], "kraken2_report")
                command_kraken2 = f"kraken2 --paired {kraken2_params_str} --report {kraken2_report} --output /dev/null --unclassified-out {kraken2_out1},{kraken2_out2} {pair[0]} {pair[1]}"
                commands.append(command_kraken2)
                bowtie2_out1 = get_bowtie2_output_filename(pair[0], "1")
                bowtie2_out2 = get_bowtie2_output_filename(pair[1], "2")
                bowtie2_report = get_report_filename(pair[0], "bowtie2_report")
                command_bowtie2 = f"bowtie2 -x {human_index} -1 {kraken2_out1} -2 {kraken2_out2} --un-conc {bowtie2_out1},{bowtie2_out2} {align_params_str} > {bowtie2_report}"
                commands.append(command_bowtie2)

        for file in single_files:
            kraken2_out = get_kraken2_output_filename(file)
            kraken2_report = get_report_filename(file, "kraken2_report")
            command_kraken2 = f"kraken2 {kraken2_params_str} --report {kraken2_report} --output /dev/null --unclassified-out {kraken2_out} {file}"
            commands.append(command_kraken2)
            bowtie2_out = get_bowtie2_output_filename(file, "")
            bowtie2_report = get_report_filename(file, "bowtie2_report")
            command bowtie2 = f"bowtie2 -x {human_index} -U {kraken2_out} --un {bowtie2_out} {align_params_str} > {bowtie2_report}"
            commands.append(command_bowtie2)

    else:
        if forward and reverse:
            kraken2_out1 = get_kraken2_output_filename(forward)
            kraken2_out2 = get_kraken2_output_filename(reverse)
            kraken2_report = get_report_filename(forward, "kraken2_report")
            command_kraken2 = f"kraken2 --paired {kraken2_params_str} --report {kraken2_report} --output /dev/null --unclassified-out {kraken2_out1},{kraken2_out2} {forward} {reverse}"
            commands.append(command_kraken2)
            bowtie2_out1 = get_bowtie2_output_filename(forward, "1")
            bowtie2_out2 = get_bowtie2_output_filename(reverse, "2")
            bowtie2_report = get_report_filename(forward, "bowtie2_report")
            command bowtie2 = f"bowtie2 -x {human_index} -1 {kraken2_out1} -2 {kraken2_out2} --un-conc {bowtie2_out1},{bowtie2_out2} {align_params_str} > {bowtie2_report}"
            commands.append(command_bowtie2)
        elif single:
            kraken2_out = get_kraken2_output_filename(single)
            kraken2_report = get_report_filename(single, "kraken2_report")
            command_kraken2 = f"kraken2 {kraken2_params_str} --report {kraken2_report} --output /dev/null --unclassified-out {kraken2_out} {single}"
            commands.append(command_kraken2)
            bowtie2_out = get_bowtie2_output_filename(single, "")
            bowtie2_report = get_report_filename(single, "bowtie2_report")
            command bowtie2 = f"bowtie2 -x {human_index} -U {kraken2_out} --un {bowtie2_out} {align_params_str} > {bowtie2_report}"
            commands.append(command bowtie2)
        else:
            print("Error: Please provide either --forward and --reverse for paired-end data, --single for single-end data, or --input-dir for a directory of FASTQ files.")
            return

    run_commands(commands, module_dir, save_intermediates)

    # Run FastQC after removal
    cleaned_files = glob.glob(os.path.join(module_dir, '*_human_bowtie2_unaligned*.fastq'))
    run_fastqc(cleaned_files, module_dir)

    # Generate MultiQC report
    generate_multiqc_report(module_dir)

def remove_host_contaminants(args):
    """
    Remove host contaminant sequences using Kraken2 and Bowtie2.
    """
    config = read_config()

    if args.host == 'custom':
        if not args.db_path_bowtie2 or not args.db_path_kraken2:
            print("Error: Please provide the path to the custom databases using --db-path-bowtie2 and --db-path-kraken2.")
            return
        bowtie2_index = args.db_path_bowtie2
        kraken2_db = args.db_path_kraken2
    else:
        if args.host == 'dog':
            kraken2_db = config['databases']['kraken2_dog_db']
            bowtie2_index = config['databases']['bowtie2_dog_index']
        elif args.host == 'cat':
            kraken2_db = config['databases']['kraken2_cat_db']
            bowtie2_index = config['databases']['bowtie2_cat_index']
        elif args.host == 'rat':
            kraken2_db = config['databases']['kraken2_rat_db']
            bowtie2_index = config['databases']['bowtie2_rat_index']
        elif args.host == 'mouse':
            kraken2_db = config['databases']['kraken2_mouse_db']
            bowtie2_index = config['databases']['bowtie2_mouse_index']
        elif args.host == 'cow':
            kraken2_db = config['databases']['kraken2_cow_db']
            bowtie2_index = config['databases']['bowtie2_cow_index']
        elif args.host == 'pig':
            kraken2_db = config['databases']['kraken2_pig_db']
            bowtie2_index = config['databases']['bowtie2_pig_index']
        elif args.host == 'horse':
            kraken2_db = config['databases']['kraken2_horse_db']
            bowtie2_index = config['databases']['bowtie2_horse_index']
        elif args.host == 'zebrafish':
            kraken2_db = config['databases']['kraken2_zebrafish_db']
            bowtie2_index = config['databases']['bowtie2_zebrafish_index']
        elif args.host == 'yeast':
            kraken2_db = config['databases']['kraken2_yeast_db']
            bowtie2_index = config['databases']['bowtie2_yeast_index']
        else:
            print("Error: Unsupported host organism.")
            return

    if not check_tool_installed('kraken2'):
        suggest_installation('kraken2')
        return
    if not check_tool_installed('bowtie2'):
        suggest_installation('bowtie2')
        return

    forward = args.forward
    reverse = args.reverse
    single = args.single
    input_dir = args.input_dir
    output_dir = args.output
    save_intermediates = args.save_intermediates
    threads = args.threads if args.threads else 1
    align_params = ["--very-sensitive-local"]  # Default parameter
    kraken2_params = []

    # Kraken2 parameters
    kraken2_params.append(f("--threads {threads}")
    if args.confidence:
        kraken2_params.append(f("--confidence {args.confidence}")
    if args.minimum_hit_groups:
        kraken2_params.append(f("--minimum-hit-groups {args.minimum_hit_groups}")
    if args.report:
        kraken2_params.append(f("--report {args.report}")
    if args.use_names:
        kraken2_params.append("--use-names")
    if args.minimum_base_quality:
        kraken2_params.append(f("--minimum-base-quality {args.minimum_base_quality}")
    if args.report_minimizer_data:
        kraken2_params.append("--report-minimizer-data")

    kraken2_params_str = " ".join(kraken2_params)

    # Bowtie2 parameters
    if args.N:
        align_params.append(f("-N {args.N}")
    if args.L:
        align_params.append(f("-L {args.L}")
    if args.i:
        align_params.append(f("-i {args.i}")
    if args.rdg:
        align_params.append(f("--rdg {args.rdg}")
    if args.rfg:
        align_params.append(f("--rfg {args.rfg}")
    if args.score_min:
        align_params.append(f("--score-min {args.score_min}")
    if args.very_fast:
        align_params.append("--very-fast")
    if args.fast:
        align_params.append("--fast")
    if args.sensitive:
        align_params.append("--sensitive")
    if args.very_sensitive:
        align_params.append("--very-sensitive")
    if args.al:
        align_params.append(f("--al {args.al}")
    if args.al_conc:
        align_params.append(f("--al-conc {args.al_conc}")
    if args.un:
        align_params.append(f("--un {args.un}")
    if args.un_conc:
        align_params.append(f("--un-conc {args.un_conc}")
    if args.k:
        align_params.append(f("-k {args.k}")
    if args.a:
        align_params.append("-a")
    if args.maxins:
        align_params.append(f("--maxins {args.maxins}")

    align_params_str = " ".join(align_params)

    project_dir = create_project_dir()
    module_dir = create_module_dir(project_dir, f'{args.host}_removal')

    commands = []

    def get_kraken2_output_filename(input_file):
        base_name = os.path.basename(input_file).replace('.fastq', '').replace('.gz', '')
        return os.path.join(module_dir, f"{base_name}_kraken2_unclassified.fastq")

    def get_bowtie2_output_filename(input_file, suffix):
        base_name = os.path.basename(input_file).replace('.fastq', '').replace('.gz', '').replace('_R1', '').replace('_R2', '').replace('_1', '').replace('_2', '')
        return os.path.join(module_dir, f"{base_name}_{args.host}_bowtie2_unaligned_{suffix}.fastq")

    def get_report_filename(input_file, suffix):
        base_name = os.path.basename(input_file).replace('.fastq', '').replace('.gz', '')
        return os.path.join(module_dir, f"{base_name}_{suffix}.txt")

    if input_dir:
        paired_files = defaultdict(list)
        single_files = []

        files = glob.glob(os.path.join(input_dir, '*.fastq*'))
        for file in files:
            if '_R1' in file or '_1' in file:
                paired_files[file.replace('_R1', '').replace('_1', '')).append(file)
            elif '_R2' in file or '_2' in file:
                paired_files[file.replace('_R2', '').replace('_2', '')).append(file)
            else:
                single_files.append(file)

        for base, pair in paired_files.items():
            if len(pair) == 2:
                kraken2_out1 = get_kraken2_output_filename(pair[0])
                kraken2_out2 = get_kraken2_output_filename(pair[1])
                kraken2_report = get_report_filename(pair[0], "kraken2_report")
                command_kraken2 = f"kraken2 --paired {kraken2_params_str} --report {kraken2_report} --output /dev/null --unclassified-out {kraken2_out1},{kraken2_out2} {pair[0]} {pair[1]}"
                commands.append(command_kraken2)
                bowtie2_out1 = get_bowtie2_output_filename(pair[0], "1")
                bowtie2_out2 = get_bowtie2_output_filename(pair[1], "2")
                bowtie2_report = get_report_filename(pair[0], "bowtie2_report")
                command_bowtie2 = f"bowtie2 -x {bowtie2_index} -1 {kraken2_out1} -2 {kraken2_out2} --un-conc {bowtie2_out1},{bowtie2_out2} {align_params_str} > {bowtie2_report}"
                commands.append(command_bowtie2)

        for file in single_files:
            kraken2_out = get_kraken2_output_filename(file)
            kraken2_report = get_report_filename(file, "kraken2_report")
            command_kraken2 = f"kraken2 {kraken2_params_str} --report {kraken2_report} --output /dev/null --unclassified-out {kraken2_out} {file}"
            commands.append(command_kraken2)
            bowtie2_out = get_bowtie2_output_filename(file, "")
            bowtie2_report = get_report_filename(file, "bowtie2_report")
            command bowtie2 = f"bowtie2 -x {bowtie2_index} -U {kraken2_out} --un {bowtie2_out} {align_params_str} > {bowtie2_report}"
            commands.append(command bowtie2)

    else:
        if forward and reverse:
            kraken2_out1 = get_kraken2_output_filename(forward)
            kraken2_out2 = get_kraken2_output_filename(reverse)
            kraken2_report = get_report_filename(forward, "kraken2_report")
            command_kraken2 = f"kraken2 --paired {kraken2_params_str} --report {kraken2_report} --output /dev/null --unclassified-out {kraken2_out1},{kraken2_out2} {forward} {reverse}"
            commands.append(command_kraken2)
            bowtie2_out1 = get_bowtie2_output_filename(forward, "1")
            bowtie2_out2 = get_bowtie2_output_filename(reverse, "2")
            bowtie2_report = get_report_filename(forward, "bowtie2_report")
            command bowtie2 = f"bowtie2 -x {bowtie2_index} -1 {kraken2_out1} -2 {kraken2_out2} --un-conc {bowtie2_out1},{bowtie2_out2} {align_params_str} > {bowtie2_report}"
            commands.append(command bowtie2)
        elif single:
            kraken2_out = get_kraken2_output_filename(single)
            kraken2_report = get_report_filename(single, "kraken2_report")
            command_kraken2 = f"kraken2 {kraken2_params_str} --report {kraken2_report} --output /dev/null --unclassified-out {kraken2_out} {single}"
            commands.append(command_kraken2)
            bowtie2_out = get_bowtie2_output_filename(single, "")
            bowtie2_report = get_report_filename(single, "bowtie2_report")
            command bowtie2 = f"bowtie2 -x {bowtie2_index} -U {kraken2_out} --un {bowtie2_out} {align_params_str} > {bowtie2_report}"
            commands.append(command bowtie2)
        else:
            print("Error: Please provide either --forward and --reverse for paired-end data, --single for single-end data, or --input-dir for a directory of FASTQ files.")
            return

    run_commands(commands, module_dir, save_intermediates)

    # Run FastQC after removal
    cleaned_files = glob.glob(os.path.join(module_dir, f'*_{args.host}_bowtie2_unaligned*.fastq'))
    run_fastqc(cleaned_files, module_dir)

    # Generate MultiQC report
    generate_multiqc_report(module_dir)

def generate_multiqc_report(output_dir):
    """
    Generate a MultiQC report from the results in the output directory.
    """
    if not check_tool_installed('multiqc'):
        suggest_installation('multiqc')
        return

    command = f"multiqc {output_dir}"
    try:
        subprocess.run(command, shell=True, check=True)
        print(f"MultiQC report generated in {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while generating the MultiQC report: {e}")

if __name__ == "__main__":
    main()
EOF

echo "Setup script has successfully created the directory structure and required files for Metagenome Cleaner."
