import os
import sys

# [Updated] Added suffixes arguments to handle specific VCF files (High-Conf only)
def merge_filtered_vcfs(input_directory, output_directory, target_sample=None, in_suffix='_filtered_1.vcf', out_suffix='_filtered_1_sum.vcf'):
    # List all files matching the input suffix
    all_files = [file for file in os.listdir(input_directory) if file.endswith(in_suffix)]

    input_files = []
    
    # Logic to handle specific sample target
    if target_sample:
        # Construct expected filename
        target_file = f"{target_sample}{in_suffix}"
        if target_file in all_files:
            input_files = [target_file]
            print(f"Processing specific sample: {target_file}")
        else:
            # Try to find file starting with sample name if exact match fails
            candidates = [f for f in all_files if f.startswith(target_sample)]
            if candidates:
                input_files = [candidates[0]]
                print(f"Processing found candidate: {candidates[0]}")
            else:
                print(f"Warning: Target file {target_file} not found in {input_directory}")
                return
    else:
        # Process all files (Batch mode)
        input_files = all_files
        print(f"Processing all {len(input_files)} files in directory.")

    for input_file in input_files:
        input_path = os.path.join(input_directory, input_file)
        output_file = input_file.replace(in_suffix, out_suffix)
        output_path = os.path.join(output_directory, output_file)

        if os.path.isfile(input_path):
            try:
                with open(input_path, 'rt') as f, open(output_path, 'w') as fsum:
                    datack = False
                    olditem = []
                    pos = -1

                    for rows in f:
                        newitem = rows.split('\t')

                        if not datack:
                            if rows[:6] == "#CHROM":
                                datack = True
                            fsum.write(rows)
                        elif datack and (len(newitem) >= 5): # Safety check for line length
                            # Logic to merge sequential SNPs (pos, pos+1)
                            if (pos + 1 == int(newitem[1])) and (len(newitem[3]) == 1) and (len(newitem[4]) == 1):
                                pos = pos + 1
                                olditem[3] = olditem[3] + newitem[3] # Merge REF
                                olditem[4] = olditem[4] + newitem[4] # Merge ALT
                            elif len(olditem) > 0:
                                stem = "\t".join(olditem)
                                fsum.write(stem)
                                pos = int(newitem[1])
                                olditem = newitem[:]
                            else:
                                pos = int(newitem[1])
                                olditem = newitem[:]
                        
                    # Write the last item
                    if len(olditem) > 0:
                        stem = "\t".join(olditem)
                        fsum.write(stem)
                        
            except Exception as e:
                print(f"Error processing {input_file}: {str(e)}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python merge_vcf.py <input_dir> <output_dir> [sample_name] [in_suffix] [out_suffix]")
        sys.exit(1)

    input_directory = sys.argv[1]
    output_directory = sys.argv[2]
    
    target_sample = None
    in_suffix = '_filtered_1.vcf'
    out_suffix = '_filtered_1_sum.vcf'

    if len(sys.argv) > 3:
        target_sample = sys.argv[3]
    
    if len(sys.argv) > 5:
        in_suffix = sys.argv[4]
        out_suffix = sys.argv[5]

    merge_filtered_vcfs(input_directory, output_directory, target_sample, in_suffix, out_suffix)