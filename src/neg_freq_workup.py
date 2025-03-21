import glob
import re
import pandas as pd


def neg_freq_workup(path_to_out_files: str, destination_path: str) -> None:
    """
    Finds all .out files recursively in a given directory and checks if they have vibrational frequency
    calculations. If so, it checks to see if any negative frequencies wer found and compiles them into
    an Excel spreadsheet.
    """
    out_file_paths = glob.glob(f"{path_to_out_files}/**/*.out", recursive=True)
    calc_name_to_neg_freqs = {}

    for out_file_path in out_file_paths:
        out_file_path = out_file_path.replace("\\", "/")

        calc_name = out_file_path.split("/")[-1][:-4]  # remove .out

        with open(out_file_path, 'r') as out_file:
            outfile_contents = out_file.read()

        # Find vibrational frequencies section, skip if it happens not found:
        try:
            regex = re.compile(r'VIBRATIONAL\ FREQUENCIES.*(\ \ \ \ \ 0:.*)NORMAL\ MODES', flags=re.DOTALL)
            search_result = regex.search(outfile_contents).group(1)
        except AttributeError:
            continue
        
        neg_freqs = ''
        split_result = search_result.splitlines('\n')
        for i in range(len(split_result)):
            try:
                freq = float(split_result[i][7:][:-7].replace(' ', ''))
                if freq < 0:
                    neg_freqs += f'  {freq}'
            except ValueError:
                continue
        
        if len(neg_freqs) == 0:
            neg_freqs = "no negative frequencies found"
        else:
            neg_freqs = neg_freqs[2:]  # remove initial spaces

        calc_name_to_neg_freqs[calc_name] = neg_freqs
    
    # Make Excel spreadsheet of the negative frequencies:
    neg_freq_df = pd.Series(calc_name_to_neg_freqs).to_frame()
    neg_freq_df = neg_freq_df.reset_index()
    neg_freq_df = neg_freq_df.sort_values(by="index")  # sort by molecule label
    neg_freq_df.to_excel(f"{destination_path}/negative_freqs.xlsx", index=False)



neg_freq_workup(
    path_to_out_files="src/data/all_bigjob/",
    destination_path="src/data/all_bigjob/",
)