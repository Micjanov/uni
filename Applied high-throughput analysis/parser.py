from argparse import ArgumentParser
from glob import glob
from typing import List, Tuple

from pandas import DataFrame, read_html, Series


def find_html_files(path) -> List:
    file_list = glob(path + "/*fastqc.html")
    return file_list


def extract_adapters(files) -> Tuple[Series, List]:
    all_adapters = DataFrame()
    for entry in files:
        tables = read_html(entry)
        adapter_sequences = tables[1].Sequence
        all_adapters = all_adapters.append(adapter_sequences)
    processed_adapters = preprocess_dataframe(all_adapters)
    stats = [
        processed_adapters.str.len().mean(),
        processed_adapters.str.len().std(),
    ]
    return processed_adapters, stats


def preprocess_dataframe(adapters) -> Series:
    na_free_adapters = adapters.dropna(axis=1)
    stacked_adapters = na_free_adapters.stack()
    duplicate_free_adapters = stacked_adapters.drop_duplicates()
    return duplicate_free_adapters


def save_to_file(filename, adapters) -> None:
    with open(filename, "w") as f:
        for index, value in adapters.iteritems():
            fasta_entry = f">{index}\n{value}\n"
            f.write(fasta_entry)


def parse_arguments():
    parser = ArgumentParser()
    parser.add_argument("input", help="directory containing the fastqc reports")
    parser.add_argument("output", help="file where to export the sequences")
    return parser.parse_args()


def main():
    args = parse_arguments()
    file_list = find_html_files(args.input)
    adapters, stats = extract_adapters(file_list)
    save_to_file(args.output, adapters)
    print(
        f"Mean of sequence length: {stats[0]}, standard deviation of sequence length {stats[1]}"
    )


if __name__ == "__main__":
    main()
