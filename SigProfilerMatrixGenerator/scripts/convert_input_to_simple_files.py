#!/usr/bin/env python3

# Author: Erik Bergstrom

# Contact: ebergstr@eng.ucsd.edu

from __future__ import print_function

import os

from SigProfilerMatrixGenerator.scripts import MutationMatrixGenerator as spm

from gzip import open as open_gz
import fcntl

def convertVCF(project, vcf_path, genome, output_path, ncbi_chrom, log_file):
    """
    Converts input vcf files into a single simple text format.

    Parameters:
             project  -> unique name given to the current samples
            vcf_path  -> path to the input vcf files
              genome  -> reference genome
     output_path  -> path to the temporary folder

    Returns:
                     snv  -> Boolean that informs whether there are SNVs present in the
                                     input files
               indel  -> Boolean that informs whether there are INDELs present
                                     in the input files

    Ouput:
            Saves a single text file to the temporary file folder

    """
    # Collect all input file names and instantiate flags
    transcript_path = (
        str(spm.reference_paths(genome)[1])
        + "/references/chromosomes/transcripts/"
        + genome
        + "/"
    )
    out_chroms = [
        x.replace("_transcripts.txt", "")
        for x in os.listdir(transcript_path)
        if not x.startswith(".")
    ]
    files = os.listdir(vcf_path)
    first_indel = True
    first_SNV = True
    snv = False
    indel = False
    out = open(log_file, "a")
    prev_line = None
    skipped_count = 0
    samples = []

    # Iterates through each file
    for file in files:
        # Skip hidden files
        if file[0] == ".":
            continue
        file_name = file.split(".")
        sample = file_name[0]
        if sample not in samples:
            samples.append(sample)
        with open(vcf_path + file) as f:
            for lines in f:
                # Skips any header lines
                if lines[0] == "#":
                    continue
                else:
                    try:
                        line = lines.strip().split()
                        if len(line) == 0:
                            continue
                        chrom = line[0]
                        if len(chrom) > 2 and genome.lower() != "ebv":
                            chrom = chrom[3:]
                        if chrom in ncbi_chrom:
                            chrom = ncbi_chrom[chrom]
                        if chrom.upper() == "M" or chrom == "mt":
                            chrom = "MT"
                        start = line[1]
                        ref = line[3]
                        mut = line[4]
                        int(start)

                    except:
                        print(
                            "The given input files do not appear to be in the correct vcf format. Skipping this file: ",
                            file,
                        )
                        break

                    # Saves SNV mutations into an SNV simple text file
                    if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
                        snv = True

                        if ref not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref == mut:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue
                        # Open all SNV files to be written to
                        if first_SNV:
                            if not os.path.exists(output_path + "SNV/"):
                                os.mkdir(output_path + "SNV/")
                            chrom_names = [
                                str(output_path)
                                + "SNV/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_files = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_names
                            ]
                            outFiles = dict(zip(out_chroms, chrom_files))
                            first_SNV = False

                        if chrom in outFiles:
                            print(
                                "\t".join([sample, chrom, start, ref, mut]),
                                file=outFiles[chrom],
                            )
                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    elif (
                        len(ref) == 2
                        and len(mut) == 2
                        and "-" not in ref
                        and "-" not in mut
                    ):
                        ref_1 = ref[0]
                        ref_2 = ref[1]
                        mut_1 = mut[0]
                        mut_2 = mut[1]
                        snv = True
                        # Check first base combination
                        if ref_1 not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut_1 not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref_1 == mut_1:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue
                        # Check second base combination
                        if ref_2 not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut_2 not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref_2 == mut_2:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if first_SNV:
                            if not os.path.exists(output_path + "SNV/"):
                                os.mkdir(output_path + "SNV/")
                            chrom_names = [
                                str(output_path)
                                + "SNV/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_files = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_names
                            ]
                            outFiles = dict(zip(out_chroms, chrom_files))
                            first_SNV = False

                        if chrom in outFiles:
                            print(
                                "\t".join([sample, chrom, start, ref_1, mut_1]),
                                file=outFiles[chrom],
                            )
                            print(
                                "\t".join(
                                    [sample, chrom, str(int(start) + 1), ref_2, mut_2]
                                ),
                                file=outFiles[chrom],
                            )
                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    # Saves INDEL mutations into an INDEL simple text file
                    else:
                        indel = True
                        if first_indel:
                            if not os.path.exists(output_path + "INDEL/"):
                                os.mkdir(output_path + "INDEL/")

                            chrom_namesI = [
                                str(output_path)
                                + "INDEL/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_filesI = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_namesI
                            ]
                            outFilesI = dict(zip(out_chroms, chrom_filesI))
                            first_indel = False

                        if chrom in outFilesI:
                            print(
                                "\t".join([sample, chrom, start, ref, mut]),
                                file=outFilesI[chrom],
                            )
                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    prev_line = line

    # Closes the output files and returns the boolean flags
    if snv:
        for files in outFiles.values():
            files.close()
    if indel:
        for files in outFilesI.values():
            files.close()
    out.close()
    return (snv, indel, skipped_count, samples)


def convertTxt(project, vcf_path, genome, output_path, ncbi_chrom, log_file):
    """
    Converts input text files into a single simple text format.

    Parameters:
             project  -> unique name given to the current samples
            vcf_path  -> path to the input text files
              genome  -> reference genome
     output_path  -> path to the temporary folder

    Returns:
                     snv  -> Boolean that informs whether there are SNVs present in the
                                     input files
               indel  -> Boolean that informs whether there are INDELs present
                                     in the input files

    Ouput:
            Saves a single text file to the temporary file folder

    """

    # Collect all input file names and instantiate flags
    transcript_path = (
        str(spm.reference_paths(genome)[1])
        + "/references/chromosomes/transcripts/"
        + genome
        + "/"
    )
    out_chroms = [
        x.replace("_transcripts.txt", "")
        for x in os.listdir(transcript_path)
        if not x.startswith(".")
    ]
    out = open(log_file, "a")
    files = os.listdir(vcf_path)
    first_indel = True
    first_SNV = True
    snv = False
    indel = False
    prev_line = None
    skipped_count = 0
    samples = []

    # Iterates through each file
    for file in files:
        if file[0] == ".":
            continue
        with open(vcf_path + file) as f:
            next(f)
            for lines in f:
                try:
                    line = lines.strip().split()
                    if len(line) == 0:
                        continue
                    sample = line[1]
                    if sample not in samples:
                        samples.append(sample)
                    chrom = line[5]
                    if len(chrom) > 2:
                        chrom = chrom[3:]
                    if chrom in ncbi_chrom:
                        chrom = ncbi_chrom[chrom]
                    if chrom.upper() == "M" or chrom == "mt":
                        chrom = "MT"
                    start = line[6]
                    end = line[7]
                    ref = line[8]
                    mut = line[9]
                    int(start)
                    int(end)

                except:
                    print(
                        "The given input files do not appear to be in the correct simple text format. Skipping this file: ",
                        file,
                    )
                    break

                # Saves SNV mutations into an SNV simple text file
                if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
                    snv = True
                    if ref not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref == mut:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if first_SNV:
                        if not os.path.exists(output_path + "SNV/"):
                            os.mkdir(output_path + "SNV/")
                        chrom_names = [
                            str(output_path)
                            + "SNV/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_files = [
                            open(file_name, "w", 10000000) for file_name in chrom_names
                        ]
                        outFiles = dict(zip(out_chroms, chrom_files))
                        first_SNV = False

                    if chrom in outFiles:
                        print(
                            "\t".join([sample, chrom, start, ref, mut]),
                            file=outFiles[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()
                elif (
                    len(ref) == 2
                    and len(mut) == 2
                    and "-" not in ref
                    and "-" not in mut
                ):
                    ref_1 = ref[0]
                    ref_2 = ref[1]
                    mut_1 = mut[0]
                    mut_2 = mut[1]
                    snv = True
                    # Check first base combination
                    if ref_1 not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut_1 not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref_1 == mut_1:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue
                    # Check second base combination
                    if ref_2 not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut_2 not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref_2 == mut_2:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if first_SNV:
                        if not os.path.exists(output_path + "SNV/"):
                            os.mkdir(output_path + "SNV/")
                        chrom_names = [
                            str(output_path)
                            + "SNV/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_files = [
                            open(file_name, "w", 10000000) for file_name in chrom_names
                        ]
                        outFiles = dict(zip(out_chroms, chrom_files))
                        first_SNV = False

                    if chrom in outFiles:
                        print(
                            "\t".join([sample, chrom, start, ref_1, mut_1]),
                            file=outFiles[chrom],
                        )
                        print(
                            "\t".join(
                                [sample, chrom, str(int(start) + 1), ref_2, mut_2]
                            ),
                            file=outFiles[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()

                # Saves INDEL mutations into an INDEL simple text file
                else:
                    indel = True
                    if first_indel:
                        if not os.path.exists(output_path + "INDEL/"):
                            os.mkdir(output_path + "INDEL/")
                        chrom_namesI = [
                            str(output_path)
                            + "INDEL/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_filesI = [
                            open(file_name, "w", 10000000) for file_name in chrom_namesI
                        ]
                        outFilesI = dict(zip(out_chroms, chrom_filesI))
                        first_indel = False

                    if chrom in outFilesI:
                        print(
                            "\t".join([sample, chrom, start, ref, mut]),
                            file=outFilesI[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()
                prev_line = line

    # Closes the output files and returns the boolean flags
    if snv:
        for files in outFiles.values():
            files.close()
        # out_snv.close()
    if indel:
        for files in outFilesI.values():
            files.close()
        # out_indel.close()
    out.close()
    return (snv, indel, skipped_count, samples)


def convertMAF(project, vcf_path, genome, output_path, ncbi_chrom, log_file):
    """
    Converts input MAF files into a single simple text format.

    Parameters:
             project  -> unique name given to the current samples
            vcf_path  -> path to the input MAF files
            genome  -> reference genome
     output_path  -> path to the temporary folder

    Returns:
                     snv  -> Boolean that informs whether there are SNVs present in the
                                     input files
               indel  -> Boolean that informs whether there are INDELs present
                                     in the input files

    Ouput:
            Saves a single text file to the temporary file folder

    """

    # Collect all input file names and instantiate flags
    out = open(log_file, "a")
    files = os.listdir(vcf_path)
    first_indel = True
    first_SNV = True
    snv = False
    indel = False
    prev_line = None
    skipped_count = 0
    samples = []

    # Iterates through each file
    transcript_path = (
        str(spm.reference_paths(genome)[1])
        + "/references/chromosomes/transcripts/"
        + genome
        + "/"
    )
    out_chroms = [
        x.replace("_transcripts.txt", "")
        for x in os.listdir(transcript_path)
        if not x.startswith(".")
    ]
    for file in files:
        header = True
        if file[0] == ".":
            continue
        name = file.split(".")
        with open(vcf_path + file) as f:
            for lines in f:
                if lines[0] == "#":
                    continue
                elif header:
                    header = False
                    continue
                try:
                    line = lines.strip().split("\t")
                    if len(line) == 0:
                        continue
                    chrom = line[4]
                    if len(chrom) > 2:
                        chrom = chrom[3:]
                    if chrom in ncbi_chrom:
                        chrom = ncbi_chrom[chrom]
                    if chrom.upper() == "M" or chrom == "mt":
                        chrom = "MT"
                    start = line[5]
                    end = line[6]
                    ref = line[10]
                    mut = line[12]
                    sample = line[15]
                    if sample not in samples:
                        samples.append(sample)
                    int(start)
                    int(end)

                except:
                    print(
                        "The given input files do not appear to be in the correct MAF format. Skipping this file: ",
                        file,
                    )
                    break

                # Saves SNV mutations into an SNV simple text file
                if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
                    snv = True
                    if ref not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref == mut:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if first_SNV:
                        if not os.path.exists(output_path + "SNV/"):
                            os.mkdir(output_path + "SNV/")
                        chrom_names = [
                            str(output_path)
                            + "SNV/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_files = [
                            open(file_name, "w", 10000000) for file_name in chrom_names
                        ]
                        outFiles = dict(zip(out_chroms, chrom_files))
                        first_SNV = False

                    if chrom in outFiles:
                        print(
                            "\t".join([sample, chrom, start, ref, mut]),
                            file=outFiles[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()

                elif (
                    len(ref) == 2
                    and len(mut) == 2
                    and "-" not in ref
                    and "-" not in mut
                ):
                    ref_1 = ref[0]
                    ref_2 = ref[1]
                    mut_1 = mut[0]
                    mut_2 = mut[1]
                    snv = True
                    # Check first base combination
                    if ref_1 not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut_1 not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref_1 == mut_1:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue
                    # Check second base combination
                    if ref_2 not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut_2 not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref_2 == mut_2:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if first_SNV:
                        if not os.path.exists(output_path + "SNV/"):
                            os.mkdir(output_path + "SNV/")
                        chrom_names = [
                            str(output_path)
                            + "SNV/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_files = [
                            open(file_name, "w", 10000000) for file_name in chrom_names
                        ]
                        outFiles = dict(zip(out_chroms, chrom_files))
                        first_SNV = False

                    if chrom in outFiles:
                        print(
                            "\t".join([sample, chrom, start, ref_1, mut_1]),
                            file=outFiles[chrom],
                        )
                        print(
                            "\t".join(
                                [sample, chrom, str(int(start) + 1), ref_2, mut_2]
                            ),
                            file=outFiles[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()

                # Saves INDEL mutations into an INDEL simple text file
                else:
                    start = str(int(line[5]) - 1)
                    indel = True
                    if first_indel:
                        if not os.path.exists(output_path + "INDEL/"):
                            os.mkdir(output_path + "INDEL/")
                        chrom_namesI = [
                            str(output_path)
                            + "INDEL/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_filesI = [
                            open(file_name, "w", 10000000) for file_name in chrom_namesI
                        ]
                        outFilesI = dict(zip(out_chroms, chrom_filesI))
                        first_indel = False

                    if chrom in outFilesI:
                        print(
                            "\t".join([sample, chrom, start, ref, mut]),
                            file=outFilesI[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()

                prev_line = line

    # Closes the output files and returns the boolean flags
    if snv:
        for files in outFiles.values():
            files.close()
        # out_snv.close()
    if indel:
        for files in outFilesI.values():
            files.close()
        # out_indel.close()
    out.close()
    return (snv, indel, skipped_count, samples)


def convertICGC(project, vcf_path, genome, output_path, ncbi_chrom, log_file):
    """
    Converts input ICGC files into a single simple text format.

    Parameters:
             project  -> unique name given to the current samples
            vcf_path  -> path to the input ICGC files
     output_path  -> path to the temporary folder
          genome  -> reference genome

    Returns:
                     snv  -> Boolean that informs whether there are SNVs present in the
                                     input files
               indel  -> Boolean that informs whether there are INDELs present
                                     in the input files

    Ouput:
            Saves a single text file to the temporary file folder

    """

    # Collect all input file names and instantiate flags
    out = open(log_file, "a")
    files = os.listdir(vcf_path)
    first_indel = True
    first_SNV = True
    snv = False
    indel = False
    prev_line = None
    skipped_count = 0
    samples = []

    # Iterates through each file
    transcript_path = (
        str(spm.reference_paths(genome)[1])
        + "/references/chromosomes/transcripts/"
        + genome
        + "/"
    )
    out_chroms = [
        x.replace("_transcripts.txt", "")
        for x in os.listdir(transcript_path)
        if not x.startswith(".")
    ]
    for file in files:
        if file[0] == ".":
            continue
        with open(vcf_path + file) as f:
            # skip the header line in the file
            next(f)
            for lines in f:
                try:
                    line = lines.strip().split("\t")
                    if len(line) == 0:
                        continue
                    sample = line[1]
                    if sample not in samples:
                        samples.append(sample)
                    icgc_sample_id = line[4]
                    chrom = line[8]
                    if len(chrom) > 2:
                        chrom = chrom[3:]
                    if chrom in ncbi_chrom:
                        chrom = ncbi_chrom[chrom]
                    if chrom.upper() == "M" or chrom == "mt":
                        chrom = "MT"
                    start = line[9]
                    end = line[10]
                    ref = line[15]
                    mut = line[16]
                    if ref == "-":
                        mut = "-" + mut
                    elif mut == "-":
                        start -= str(int(start) - 1)
                        ref = "-" + ref
                    int(start)
                    int(end)
                except:
                    print(
                        "The given input files do not appear to be in the correct ICGC format."
                    )
                    break

                # Saves SNV mutations into an SNV simple text file
                if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
                    snv = True
                    if ref not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref == mut:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if first_SNV:
                        if not os.path.exists(output_path + "SNV/"):
                            os.mkdir(output_path + "SNV/")
                        chrom_names = [
                            str(output_path)
                            + "SNV/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_files = [
                            open(file_name, "w", 10000000) for file_name in chrom_names
                        ]
                        outFiles = dict(zip(out_chroms, chrom_files))
                        first_SNV = False

                    if chrom in outFiles:
                        print(
                            "\t".join([sample, chrom, start, ref, mut]),
                            file=outFiles[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()

                elif (
                    len(ref) == 2
                    and len(mut) == 2
                    and "-" not in ref
                    and "-" not in mut
                ):
                    ref_1 = ref[0]
                    ref_2 = ref[1]
                    mut_1 = mut[0]
                    mut_2 = mut[1]
                    snv = True
                    # Check first base combination
                    if ref_1 not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut_1 not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref_1 == mut_1:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue
                    # Check second base combination
                    if ref_2 not in "ACGT-":
                        print(
                            "The ref base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if mut_2 not in "ACGT-":
                        print(
                            "The mutation base is not recognized. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if ref_2 == mut_2:
                        print(
                            "The ref base appears to match the mutated base. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if line == prev_line:
                        print(
                            "There appears to be a duplicate single base substitution. Skipping this mutation: "
                            + chrom
                            + " "
                            + str(start)
                            + " "
                            + ref
                            + " "
                            + mut,
                            file=out,
                        )
                        out.flush()
                        skipped_count += 1
                        continue

                    if first_SNV:
                        if not os.path.exists(output_path + "SNV/"):
                            os.mkdir(output_path + "SNV/")
                        chrom_names = [
                            str(output_path)
                            + "SNV/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_files = [
                            open(file_name, "w", 10000000) for file_name in chrom_names
                        ]
                        outFiles = dict(zip(out_chroms, chrom_files))
                        first_SNV = False

                    if chrom in outFiles:
                        print(
                            "\t".join([sample, chrom, start, ref_1, mut_1]),
                            file=outFiles[chrom],
                        )
                        print(
                            "\t".join(
                                [sample, chrom, str(int(start) + 1), ref_2, mut_2]
                            ),
                            file=outFiles[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()

                # Saves INDEL mutations into an INDEL simple text file
                else:
                    indel = True
                    if first_indel:
                        if not os.path.exists(output_path + "INDEL/"):
                            os.mkdir(output_path + "INDEL/")
                        chrom_namesI = [
                            str(output_path)
                            + "INDEL/"
                            + out_chroms[i]
                            + "_"
                            + project
                            + ".genome"
                            for i in range(0, len(out_chroms))
                        ]
                        chrom_filesI = [
                            open(file_name, "w", 10000000) for file_name in chrom_namesI
                        ]
                        outFilesI = dict(zip(out_chroms, chrom_filesI))
                        first_indel = False

                    if chrom in outFilesI:
                        print(
                            "\t".join([sample, chrom, start, ref, mut]),
                            file=outFilesI[chrom],
                        )
                    else:
                        print(
                            chrom
                            + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                            file=out,
                        )
                        out.flush()

                prev_line = line

    # Closes the output files and returns the boolean flags
    if snv:
        for files in outFiles.values():
            files.close()
        # out_snv.close()
    if indel:
        for files in outFilesI.values():
            files.close()
        # out_indel.close()
    out.close()
    return (snv, indel, skipped_count, samples)


def convertMutliSampleVCF(project, vcf_path, genome, output_path, ncbi_chrom, log_file, binary_gz=False):
    """
    Alternative function to handle vcf files with calls from multiple samples.
    Tested on vcf files from the 1000 Genomes Project.
    Converts input vcf files into a single simple text format.
    Able to read vcf.gz files using the gzip module if binary_gz is set True.

    Parameters:
            project  -> unique name given to the current samples
            vcf_path  -> path to the input vcf files
            genome  -> reference genome
            output_path  -> path to the temporary folder


    Returns:
                     snv  -> Boolean that informs whether there are SNVs present in the
                                     input files
               indel  -> Boolean that informs whether there are INDELs present
                                     in the input files

    Ouput:
            Saves a single text file to the temporary file folder

    """
    # Collect all input file names and instantiate flags
    transcript_path = (
        str(spm.reference_paths(genome)[1])
        + "/references/chromosomes/transcripts/"
        + genome
        + "/"
    )
    out_chroms = [
        x.replace("_transcripts.txt", "")
        for x in os.listdir(transcript_path)
        if not x.startswith(".")
    ]
    files = os.listdir(vcf_path)
    first_indel = True
    first_SNV = True
    snv = False
    indel = False
    out = open(log_file, "a")
    prev_line = None
    skipped_count = 0
    samples = []

    # Iterates through each file
    for file in files:
        # Skip hidden files
        if file[0] == ".":
            continue
        file_name = file.split(".")

        if binary_gz is False:
            mode = "r"
            open_cmd = open
        else:
            open_cmd = open_gz
            mode = "rt"
            
        with open_cmd(vcf_path + file, mode) as f:
            for lines in f:
                # Skips any header lines but use column names to retrieve samples
                if lines[0:2] == "##":
                    continue
                elif lines[0] == "#":
                    line = lines.strip().split()
                    file_samples = line[9:]
                    if len(file_samples) == 0:
                        print(
                            "No samples could be retrieved from the vcf header. Skipping this file: ",
                            file,
                        )
                        break
                    [samples.append(s) for s in file_samples if s not in samples]
                    continue
                else:
                    try:
                        line = lines.strip().split()
                        if len(line) == 0:
                            continue
                        chrom = line[0]
                        if len(chrom) > 2 and genome.lower() != "ebv":
                            chrom = chrom[3:]
                        if chrom in ncbi_chrom:
                            chrom = ncbi_chrom[chrom]
                        if chrom.upper() == "M" or chrom == "mt":
                            chrom = "MT"
                        start = line[1]
                        ref = line[3]
                        mut = line[4]
                        int(start)
                        # Check that genotype is defined in format column
                        form = line[8]
                        if form.split(":")[0] != "GT":
                            print(
                                "The format column value does not start with GT as expected. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut
                                + " "
                                + form,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                    except:
                        print(
                            "The given input files do not appear to be in the correct vcf format. Skipping this file: ",
                            file,
                        )
                        break

                    # Saves SNV mutations into an SNV simple text file
                    if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
                        snv = True

                        if ref not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref == mut:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue
                        # Open all SNV files to be written to
                        if first_SNV:
                            if not os.path.exists(output_path + "SNV/"):
                                os.mkdir(output_path + "SNV/")
                            chrom_names = [
                                str(output_path)
                                + "SNV/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_files = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_names
                            ]
                            outFiles = dict(zip(out_chroms, chrom_files))
                            first_SNV = False

                        if chrom in outFiles:

                            for sample, sample_info in zip(file_samples, line[9:]):
                                gt = sample_info.split(":")[0]

                                # Skip sample entries without the variant
                                if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                    continue
                                
                                print(
                                    "\t".join([sample, chrom, start, ref, mut]),
                                    file=outFiles[chrom],
                                )
                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    elif (
                        len(ref) == 2
                        and len(mut) == 2
                        and "-" not in ref
                        and "-" not in mut
                    ):
                        ref_1 = ref[0]
                        ref_2 = ref[1]
                        mut_1 = mut[0]
                        mut_2 = mut[1]
                        snv = True
                        # Check first base combination
                        if ref_1 not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut_1 not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref_1 == mut_1:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue
                        # Check second base combination
                        if ref_2 not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut_2 not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref_2 == mut_2:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if first_SNV:
                            if not os.path.exists(output_path + "SNV/"):
                                os.mkdir(output_path + "SNV/")
                            chrom_names = [
                                str(output_path)
                                + "SNV/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_files = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_names
                            ]
                            outFiles = dict(zip(out_chroms, chrom_files))
                            first_SNV = False

                        if chrom in outFiles:

                            for sample, sample_info in zip(file_samples, line[9:]):
                                gt = sample_info.split(":")[0]

                                # Skip sample entries without the variant
                                if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                    continue
                                
                                print(
                                    "\t".join([sample, chrom, start, ref_1, mut_1]),
                                    file=outFiles[chrom],
                                )
                                print(
                                    "\t".join(
                                        [sample, chrom, str(int(start) + 1), ref_2, mut_2]
                                    ),
                                    file=outFiles[chrom],
                            )
                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    # Saves INDEL mutations into an INDEL simple text file
                    else:
                        indel = True
                        if first_indel:
                            if not os.path.exists(output_path + "INDEL/"):
                                os.mkdir(output_path + "INDEL/")

                            chrom_namesI = [
                                str(output_path)
                                + "INDEL/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_filesI = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_namesI
                            ]
                            outFilesI = dict(zip(out_chroms, chrom_filesI))
                            first_indel = False

                        if chrom in outFilesI:

                            for sample, sample_info in zip(file_samples, line[9:]):
                                gt = sample_info.split(":")[0]

                                # Skip sample entries without the variant
                                if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                    continue
                                
                                print(
                                    "\t".join([sample, chrom, start, ref, mut]),
                                    file=outFilesI[chrom],
                                )

                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    prev_line = line

    # Closes the output files and returns the boolean flags
    if snv:
        for files in outFiles.values():
            files.close()
    if indel:
        for files in outFilesI.values():
            files.close()
    out.close()
    return (snv, indel, skipped_count, samples)


def _wip_convertMutliSampleVCF(project, vcf_path, genome, output_path, ncbi_chrom, log_file, binary_gz=False):
    """

    wip: multiprocessing support

    Alternative function to handle vcf files with calls from multiple samples.
    Tested on vcf files from the 1000 Genomes Project.
    Converts input vcf files into a single simple text format.
    Able to read vcf.gz files using the gzip module if binary_gz is set True.

    Parameters:
            project  -> unique name given to the current samples
            vcf_path  -> path to the input vcf files
            genome  -> reference genome
            output_path  -> path to the temporary folder


    Returns:
                     snv  -> Boolean that informs whether there are SNVs present in the
                                     input files
               indel  -> Boolean that informs whether there are INDELs present
                                     in the input files

    Ouput:
            Saves a single text file to the temporary file folder

    """
    # Collect all input file names and instantiate flags
    transcript_path = (
        str(spm.reference_paths(genome)[1])
        + "/references/chromosomes/transcripts/"
        + genome
        + "/"
    )
    out_chroms = [
        x.replace("_transcripts.txt", "")
        for x in os.listdir(transcript_path)
        if not x.startswith(".")
    ]
    files = os.listdir(vcf_path)
    first_indel = True
    first_SNV = True
    snv = False
    indel = False
    out = open(log_file, "a")
    prev_line = None
    skipped_count = 0
    samples = []

    # Define function to handle single vcf file
    # to be used by multiprocessing module

    # handle somewhere: update samples 
    #   [samples.append(s) for s in file_samples if s not in samples]

    def parse_vcf(
            file:str,
            binary_gz=binary_gz,

        ):

        # Skip hidden files
        if file[0] == ".":
            return
        file_name = file.split(".")

        if binary_gz is False:
            mode = "r"
            open_cmd = open
        else:
            open_cmd = open_gz
            mode = "rt"
        
        # Save the first and only chromosome value of this file
        file_chr = None

        with open_cmd(vcf_path + file, mode) as f:
            for lines in f:
                # Skips any header lines but use column names to retrieve samples
                if lines[0:2] == "##":
                    continue
                elif lines[0] == "#":
                    line = lines.strip().split()
                    file_samples = line[9:]
                    if len(file_samples) == 0:
                        print(
                            "No samples could be retrieved from the vcf header. Skipping this file: ",
                            file,
                        )
                        break
                    continue
                else:
                    try:
                        line = lines.strip().split()
                        if len(line) == 0:
                            continue
                        chrom = line[0]
                        if len(chrom) > 2 and genome.lower() != "ebv":
                            chrom = chrom[3:]
                        if chrom in ncbi_chrom:
                            chrom = ncbi_chrom[chrom]
                        if chrom.upper() == "M" or chrom == "mt":
                            chrom = "MT"

                        if file_chr is None:
                            file_chr = chrom

                        start = line[1]
                        ref = line[3]
                        mut = line[4]
                        int(start)
                        # Check that genotype is defined in format column
                        form = line[8]
                        if form.split(":")[0] != "GT":
                            # Lock file as to not generate corruption when writing by multiple processes
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The format column value does not start with GT as expected. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut
                                + " "
                                + form,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                    except:
                        print(
                            "The given input files do not appear to be in the correct vcf format. Skipping this file: ",
                            file,
                        )
                        break

                    # Saves SNV mutations into an SNV simple text file
                    if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
                        snv = True

                        if ref not in "ACGT-":
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if mut not in "ACGT-":
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if ref == mut:
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue
                        # Open all SNV files to be written to
                        # Open only the one required SNV file
                        if first_SNV:
                            os.makedirs(output_path + "SNV/", exist_ok=True)
                            chrom_names = [
                                str(output_path)
                                + "SNV/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_files = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_names
                            ]
                            outFiles = dict(zip(out_chroms, chrom_files))
                            first_SNV = False

                        if chrom in outFiles:

                            for sample, sample_info in zip(file_samples, line[9:]):
                                gt = sample_info.split(":")[0]

                                # Skip sample entries without the variant
                                if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                    continue
                                
                                print(
                                    "\t".join([sample, chrom, start, ref, mut]),
                                    file=outFiles[chrom],
                                )
                        else:
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)

                    elif (
                        len(ref) == 2
                        and len(mut) == 2
                        and "-" not in ref
                        and "-" not in mut
                    ):
                        ref_1 = ref[0]
                        ref_2 = ref[1]
                        mut_1 = mut[0]
                        mut_2 = mut[1]
                        snv = True
                        # Check first base combination
                        if ref_1 not in "ACGT-":
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if mut_1 not in "ACGT-":
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if ref_1 == mut_1:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue
                        # Check second base combination
                        if ref_2 not in "ACGT-":
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if mut_2 not in "ACGT-":
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if ref_2 == mut_2:
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)
                            skipped_count += 1
                            continue

                        if first_SNV:
                            if not os.path.exists(output_path + "SNV/"):
                                os.mkdir(output_path + "SNV/")
                            chrom_names = [
                                str(output_path)
                                + "SNV/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_files = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_names
                            ]
                            outFiles = dict(zip(out_chroms, chrom_files))
                            first_SNV = False

                        if chrom in outFiles:

                            for sample, sample_info in zip(file_samples, line[9:]):
                                gt = sample_info.split(":")[0]

                                # Skip sample entries without the variant
                                if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                    continue
                                
                                print(
                                    "\t".join([sample, chrom, start, ref_1, mut_1]),
                                    file=outFiles[chrom],
                                )
                                print(
                                    "\t".join(
                                        [sample, chrom, str(int(start) + 1), ref_2, mut_2]
                                    ),
                                    file=outFiles[chrom],
                            )
                        else:
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            fcntl.flock(out, fcntl.LOCK_UN)
                            out.flush()

                    # Saves INDEL mutations into an INDEL simple text file
                    else:
                        indel = True

                        # Open only the one required INDEL file
                        if first_indel:
                            os.makedirs(output_path + "INDEL/", exist_ok=True)

                            chrom_namesI = [
                                str(output_path)
                                + "INDEL/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_filesI = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_namesI
                            ]
                            outFilesI = dict(zip(out_chroms, chrom_filesI))
                            first_indel = False

                        if chrom in outFilesI:

                            for sample, sample_info in zip(file_samples, line[9:]):
                                gt = sample_info.split(":")[0]

                                # Skip sample entries without the variant
                                if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                    continue
                                
                                print(
                                    "\t".join([sample, chrom, start, ref, mut]),
                                    file=outFilesI[chrom],
                                )

                        else:
                            fcntl.flock(out, fcntl.LOCK_EX)
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()
                            fcntl.flock(out, fcntl.LOCK_UN)

                    prev_line = line

                # Close opened SNV and INDEL files
    
    return (snv, indel, skipped_count, samples)


def my_convertVCF(project, vcf_path, genome, output_path, ncbi_chrom, log_file):
    """
    
    Converts input vcf files into a single simple text format.

    MODIFIED: Checks the genotype field of each vcf entry.
    If variant is present on both alleles, enters the entry twice to the output file.
    Will skip mutations with other integers than 0 and 1 in the genotype field.
    Only set up to work with haploid/diploid inputs.

    Parameters:
             project  -> unique name given to the current samples
            vcf_path  -> path to the input vcf files
              genome  -> reference genome
     output_path  -> path to the temporary folder

    Returns:
                     snv  -> Boolean that informs whether there are SNVs present in the
                                     input files
               indel  -> Boolean that informs whether there are INDELs present
                                     in the input files

    Ouput:
            Saves a single text file to the temporary file folder

    """
    # Collect all input file names and instantiate flags
    transcript_path = (
        str(spm.reference_paths(genome)[1])
        + "/references/chromosomes/transcripts/"
        + genome
        + "/"
    )
    out_chroms = [
        x.replace("_transcripts.txt", "")
        for x in os.listdir(transcript_path)
        if not x.startswith(".")
    ]
    files = os.listdir(vcf_path)
    first_indel = True
    first_SNV = True
    snv = False
    indel = False
    out = open(log_file, "a")
    prev_line = None
    skipped_count = 0
    samples = []

    # Iterates through each file
    for file in files:
        # Skip hidden files
        if file[0] == ".":
            continue
        file_name = file.split(".")
        sample = file_name[0]
        if sample not in samples:
            samples.append(sample)
        with open(vcf_path + file) as f:
            for lines in f:
                # Skips any header lines
                if lines[0] == "#":
                    continue
                else:
                    try:
                        line = lines.strip().split()
                        if len(line) == 0:
                            continue
                        chrom = line[0]
                        if len(chrom) > 2 and genome.lower() != "ebv":
                            chrom = chrom[3:]
                        if chrom in ncbi_chrom:
                            chrom = ncbi_chrom[chrom]
                        if chrom.upper() == "M" or chrom == "mt":
                            chrom = "MT"
                        start = line[1]
                        ref = line[3]
                        mut = line[4]
                        int(start)
                        # Check that genotype is defined in format column
                        form = line[8]
                        if form.split(":")[0] != "GT":
                            print(
                                "The format column value does not start with GT as expected. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut
                                + " "
                                + form,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                    except:
                        print(
                            "The given input files do not appear to be in the correct vcf format. Skipping this file: ",
                            file,
                        )
                        break

                    # Saves SNV mutations into an SNV simple text file
                    if len(ref) == 1 and len(mut) == 1 and ref != "-" and mut != "-":
                        snv = True

                        if ref not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref == mut:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue
                        # Open all SNV files to be written to
                        if first_SNV:
                            if not os.path.exists(output_path + "SNV/"):
                                os.mkdir(output_path + "SNV/")
                            chrom_names = [
                                str(output_path)
                                + "SNV/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_files = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_names
                            ]
                            outFiles = dict(zip(out_chroms, chrom_files))
                            first_SNV = False

                        if chrom in outFiles:
                            gt = line[9].split(":")[0]
                            def print_entry(sample=sample, chrom=chrom, start=start, ref=ref, mut=mut, outFiles=outFiles):
                                print(
                                    "\t".join([sample, chrom, start, ref, mut]),
                                    file=outFiles[chrom],
                                )
                                return
                            # Skip sample entries without the variant
                            if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                print(
                                    "The vcf file contains entries for genotypes that were not found in the sample. Skipping mutation: "
                                    + chrom
                                    + " "
                                    + str(start)
                                    + " "
                                    + ref
                                    + " "
                                    + mut
                                    + " with genotype value "
                                    + gt,
                                    file=out,
                                )
                                out.flush()
                                skipped_count += 1
                                continue
                            elif str(gt) in ["1","1/0","0/1","1|0","0|1"]:
                                print_entry()
                            elif str(gt) in ["1/1","1|1"]:
                                print_entry()
                                print_entry()
                            else:
                                print(
                                    "The vcf file contains genotypes that are not supported. Skipping mutation: "
                                    + chrom
                                    + " "
                                    + str(start)
                                    + " "
                                    + ref
                                    + " "
                                    + mut
                                    + " with genotype value "
                                    + gt,
                                    file=out,
                                )                          
                                out.flush()
                                skipped_count += 1
                                continue
                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    elif (
                        len(ref) == 2
                        and len(mut) == 2
                        and "-" not in ref
                        and "-" not in mut
                    ):
                        ref_1 = ref[0]
                        ref_2 = ref[1]
                        mut_1 = mut[0]
                        mut_2 = mut[1]
                        snv = True
                        # Check first base combination
                        if ref_1 not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut_1 not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref_1 == mut_1:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue
                        # Check second base combination
                        if ref_2 not in "ACGT-":
                            print(
                                "The ref base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if mut_2 not in "ACGT-":
                            print(
                                "The mutation base is not recognized. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if ref_2 == mut_2:
                            print(
                                "The ref base appears to match the mutated base. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if line == prev_line:
                            print(
                                "There appears to be a duplicate single base substitution. Skipping this mutation: "
                                + chrom
                                + " "
                                + str(start)
                                + " "
                                + ref
                                + " "
                                + mut,
                                file=out,
                            )
                            out.flush()
                            skipped_count += 1
                            continue

                        if first_SNV:
                            if not os.path.exists(output_path + "SNV/"):
                                os.mkdir(output_path + "SNV/")
                            chrom_names = [
                                str(output_path)
                                + "SNV/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_files = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_names
                            ]
                            outFiles = dict(zip(out_chroms, chrom_files))
                            first_SNV = False

                        if chrom in outFiles:
                            gt = line[9].split(":")[0]
                            def print_entry(sample=sample, chrom=chrom, start=start, ref=ref, mut=mut, outFiles=outFiles):
                                print(
                                    "\t".join([sample, chrom, start, ref_1, mut_1]),
                                    file=outFiles[chrom],
                                )
                                print(
                                    "\t".join(
                                        [sample, chrom, str(int(start) + 1), ref_2, mut_2]
                                    ),
                                    file=outFiles[chrom],
                                )
                                return
                            # Skip sample entries without the variant
                            if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                print(
                                    "The vcf file contains entries for genotypes that were not found in the sample. Skipping mutation: "
                                    + chrom
                                    + " "
                                    + str(start)
                                    + " "
                                    + ref
                                    + " "
                                    + mut
                                    + " with genotype value "
                                    + gt,
                                    file=out,
                                )
                                out.flush()
                                skipped_count += 1
                                continue
                            elif str(gt) in ["1","1/0","0/1","1|0","0|1"]:
                                print_entry()
                            elif str(gt) in ["1/1","1|1"]:
                                print_entry()
                                print_entry()
                            else:
                                print(
                                    "The vcf file contains genotypes that are not supported. Skipping mutation: "
                                    + chrom
                                    + " "
                                    + str(start)
                                    + " "
                                    + ref
                                    + " "
                                    + mut
                                    + " with genotype value "
                                    + gt,
                                    file=out,
                                )                          
                                out.flush()
                                skipped_count += 1
                                continue
                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    # Saves INDEL mutations into an INDEL simple text file
                    else:
                        indel = True
                        if first_indel:
                            if not os.path.exists(output_path + "INDEL/"):
                                os.mkdir(output_path + "INDEL/")

                            chrom_namesI = [
                                str(output_path)
                                + "INDEL/"
                                + out_chroms[i]
                                + "_"
                                + project
                                + ".genome"
                                for i in range(0, len(out_chroms))
                            ]
                            chrom_filesI = [
                                open(file_name, "w", 10000000)
                                for file_name in chrom_namesI
                            ]
                            outFilesI = dict(zip(out_chroms, chrom_filesI))
                            first_indel = False

                        if chrom in outFilesI:
                            gt = line[9].split(":")[0]
                            def print_entry(sample=sample, chrom=chrom, start=start, ref=ref, mut=mut, outFilesI=outFilesI):
                                print(
                                    "\t".join([sample, chrom, start, ref, mut]),
                                    file=outFilesI[chrom],
                                )
                                return
                            # Skip sample entries without the variant
                            if str(gt) in [".","./.",".|.","0","0/0","0|0"]:
                                print(
                                    "The vcf file contains entries for genotypes that were not found in the sample. Skipping mutation: "
                                    + chrom
                                    + " "
                                    + str(start)
                                    + " "
                                    + ref
                                    + " "
                                    + mut
                                    + " with genotype value "
                                    + gt,
                                    file=out,
                                )
                                out.flush()
                                skipped_count += 1
                                continue
                            elif str(gt) in ["1","1/0","0/1","1|0","0|1"]:
                                print_entry()
                            elif str(gt) in ["1/1","1|1"]:
                                print_entry()
                                print_entry()
                            else:
                                print(
                                    "The vcf file contains genotypes that are not supported. Skipping mutation: "
                                    + chrom
                                    + " "
                                    + str(start)
                                    + " "
                                    + ref
                                    + " "
                                    + mut
                                    + " with genotype value "
                                    + gt,
                                    file=out,
                                )                          
                                out.flush()
                                skipped_count += 1
                                continue
                        else:
                            print(
                                chrom
                                + " is not supported. You will need to download that chromosome and create the required files. Continuing with the matrix generation...",
                                file=out,
                            )
                            out.flush()

                    prev_line = line

    # Closes the output files and returns the boolean flags
    if snv:
        for files in outFiles.values():
            files.close()
    if indel:
        for files in outFilesI.values():
            files.close()
    out.close()
    return (snv, indel, skipped_count, samples)