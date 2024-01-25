#!/bin/bash
#
# Copyright (C): 2020-2021 - Gert Hulselmans
#
# Purpose: Sort VCF file in the same order as the BAM file, so it can be used with popscle.



# Function to check if any of the programs in a pipe failed.
check_exit_codes () {
    local GET_PIPESTATUS="${PIPESTATUS[@]}";
    local exit_code;

    for exit_code in ${GET_PIPESTATUS} ; do
        if [ ${exit_code} -ne 0 ] ; then
             return ${exit_code};
        fi
    done

    return 0;
}



# Check if necessary programs are installed.
check_if_programs_exists () {
    local exit_code=0;

    # Check if awk is installed.
    if ! type awk > /dev/null 2>&1 ; then
        printf 'Error: "awk" could not be found in PATH.\n' > /dev/stderr;
        exit_code=2;
    fi

    # Check if bcftools is installed.
    if ! type bcftools > /dev/null 2>&1 ; then
        printf 'Error: "bcftools" could not be found in PATH.\n' > /dev/stderr;
        exit_code=2;
    fi

    # Check if samtools is installed.
    if ! type samtools > /dev/null 2>&1 ; then
        printf 'Error: "samtools" could not be found in PATH.\n' > /dev/stderr;
        exit_code=2;
    fi

    return ${exit_code};
}



# Get order of the contigs (chromosomes) and their length from the BAM header.
get_contig_order_from_bam () {
    local bam_input_file="${1}";
    local output_type="${2}";

    if [ ${#@} -ne 2 ] ; then
        printf 'Usage: get_contig_order_from_bam BAM_file output_type\n\n';
        printf 'Arguments:\n';
        printf '  - BAM_file: BAM file from which to get the contig order and contig lengths.\n';
        printf '  - output_type:\n';
        printf '      - "names":        Return contig names.\n';
        printf '      - "chrom_sizes":  Return contig names and contig lengths.\n';
        printf '      - "vcf":          Return VCF header section for contigs.\n\n';
        return 1;
    fi

    case "${output_type}" in
        'names')
            ;;
        'chrom_sizes')
            ;;
        'vcf')
            ;;
        *)
            printf 'Error: output_type "%s" is not supported.\n' "${output_type}" > /dev/stderr;
            return 1;
            ;;
    esac

    check_if_programs_exists || return $?;

    # Get the order of the contigs from the BAM header.
    samtools view -H "${bam_input_file}" \
      | awk \
            -F '\t' \
            -v output_type="${output_type}" \
            '
            {
                # Only look at sequence header fields.
                if ($1 == "@SQ") {
                    contig_idx += 1;
                    contig_name = "";
                    contig_length = "";

                    # Extract contig (chromosome) name and contig (chromosome) length.
                    for (i = 2; i <= NF; i++) {
                        if ($i ~ /^SN:/) {
                            contig_name = substr($i, 4);
                        }

                        if ($i ~ /^LN:/) {
                            contig_length = substr($i, 4);
                        }

                        # Create contig order to name and contig order to length and vcf contig appings.
                        contig_idx_to_name[contig_idx] = contig_name;
                        contig_idx_to_length[contig_idx] = contig_length;
                        contig_idx_to_vcf_contig[contig_idx] = sprintf("##contig=<ID=%s,length=%s>", contig_name, contig_length);
                    }
                }
            } END {
                if (contig_idx == 0) {
                    printf "Error: No \"@SQ\" header line found in BAM file.\n" > "/dev/stderr";
                    exit(1);
                } else if (output_type == "names") {
                    contig_names = "";

                    for (contig_idx = 1; contig_idx <= length(contig_idx_to_name); contig_idx++) {
                        contig_names = contig_names " " contig_idx_to_name[contig_idx];
                    }

                    # Print all contig names (without leading space).
                    print substr(contig_names, 2);
                } else if (output_type == "chrom_sizes") {
                    # Print all contig names with their length in a TAB separated fashion.
                    for (contig_idx = 1; contig_idx <= length(contig_idx_to_name); contig_idx++) {
                        print contig_idx_to_name[contig_idx] "\t" contig_idx_to_length[contig_idx];
                    }
                } else if (output_type == "vcf") {
                    # Print VCF header section for contigs.
                    for (contig_idx = 1; contig_idx <= length(contig_idx_to_vcf_contig); contig_idx++) {
                        print contig_idx_to_vcf_contig[contig_idx];
                    }
                }
            }'

      check_exit_codes;

      return $?;
}



# Sort VCF file in the same order as the BAM file, so it can be used with popscle.
sort_vcf_same_as_bam () {
    local bam_input_file="${1}";
    local vcf_input_file="${2}";
    local vcf_type="${3:-v}";

    if [ ${#@} -lt 2 ] ; then
        printf 'Usage: sort_vcf_same_as_bam BAM_file VCF_file [VCF_type]\n\n';
        printf 'Arguments:\n';
        printf '  - BAM_file: BAM file from which to get the contig order to sort the VCF file.\n';
        printf '  - VCF_file: VCF file to sort by contig order as defined in the BAM file.\n';
        printf '  - VCF_type: VCF ouput file type (default: same as input VCF file type):\n';
        printf '              v: uncompressed VCF, z: compressed VCF,\n';
        printf '              u: uncompressed BCF, b: compressed BCF\n\n';
        printf 'Purpose:\n';
        printf '  Sort VCF file in the same order as the BAM file, so it can be used with popscle.\n\n';
        return 1;
    fi

    check_if_programs_exists || return $?;

    # If VCF type is not specified, try to guess it from the filename extension.
    if [ ${#@} -eq 2 ] ; then
        if [ "${vcf_input_file%.vcf.gz}" != "${vcf_input_file}" ] ; then
            vcf_type='z';
        elif [ "${vcf_input_file%.bcf}" != "${vcf_input_file}" ] ; then
             vcf_type='b';
        fi
    fi

    # Sort VCF file by same chromosome order as BAM file.
    cat <(
          # Create new VCF header:
          #   - Get VCF header of VCF input file.
          #   - Remove all contig header lines and "#CHROM" line from the VCF header.
          #   - Append contig headers in the order they appear in the input BAM file.
          #   - Add "#CHROM" line as last line of the new VCF header.
          bcftools view -h "${vcf_input_file}" \
            | awk \
                '
                {
                    if ($1 !~ /^##contig=/ && $1 !~ /^#CHROM/) {
                        # Remove all contig header lines and "#CHROM" line.
                        print $0;
                    }
                }' \
            | cat \
                - \
                <(get_contig_order_from_bam "${bam_input_file}" 'vcf') \
                <(bcftools view -h "${vcf_input_file}" | tail -n 1) \
        ) \
        <(bcftools view -H -O v "${vcf_input_file}") \
      | bcftools sort -O "${vcf_type}";

    check_exit_codes;

    return $?;
}



sort_vcf_same_as_bam "${@}";
