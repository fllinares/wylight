from subprocess import call
import sys


def parse_lamp_output(input_filename):
    f = open(input_filename)  # Open file
    f.readline()  # Remove dummy line
    # Read from file the corrected significance threshold and testability minimum
    # support
    sig_th = f.readline().strip().split()[-1]
    lcm_th = f.readline().strip().split()[-1]
    f.close()  # Close file
    return (sig_th, lcm_th)


def execute_lamp(output_basefilename, target_fwer, input_class_labels_file, input_transactions_file,
                 lamp_path, pval_enum_path):
    output_basefilename_lamp = output_basefilename + '_lamp'
    output_basefilename_sig = output_basefilename + '_sig'
    print "Finding corrected significance threshold..."
    call([lamp_path, output_basefilename_lamp, target_fwer, input_class_labels_file, input_transactions_file])
    (sig_th, lcm_th) = parse_lamp_output(output_basefilename_lamp+'_results.txt')
    print "Finding significant itemsets..."
    call([pval_enum_path, output_basefilename_sig, sig_th, lcm_th, input_class_labels_file, input_transactions_file])
    print "Done!"


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print "Incorrect usage!"
        print "INPUT ARGUMENTS: output_basefilename target_fwer input_class_labels_file input_transctions file lamp_path pval_enum_path"
    else:
        execute_lamp(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
