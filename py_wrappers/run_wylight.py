from subprocess import call
import sys


def parse_wylight_output(input_filename):
    f = open(input_filename)  # Open file
    f.readline()  # Remove dummy line
    # Read from file the corrected significance threshold and testability minimum
    # support
    sig_th = f.readline().strip().split()[-1]
    f.readline()  # Remove dummy line
    lcm_th = f.readline().strip().split()[-1]
    f.close()  # Close file
    return (sig_th, lcm_th)


def execute_wylight(output_basefilename, n_perm, target_fwer, input_class_labels_file, input_transactions_file,
                    rand_seed, wylight_path, pval_enum_path):
    output_basefilename_wylight = output_basefilename + '_wylight'
    output_basefilename_sig = output_basefilename + '_sig'
    print "Finding corrected significance threshold..."
    call([wylight_path, output_basefilename_wylight, n_perm, target_fwer, input_class_labels_file,
          input_transactions_file, rand_seed])
    (sig_th, lcm_th) = parse_wylight_output(output_basefilename_wylight+'_results.txt')
    print "Finding significant itemsets..."
    call([pval_enum_path, output_basefilename_sig, sig_th, lcm_th, input_class_labels_file, input_transactions_file])
    print "Done!"


if __name__ == "__main__":
    if len(sys.argv) != 9:
        print "Incorrect usage!"
        print "INPUT ARGUMENTS: output_basefilename n_perm target_fwer input_class_labels_file input_transctions file rand_seed wylight_path pval_enum_path"
    else:
        execute_wylight(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])
