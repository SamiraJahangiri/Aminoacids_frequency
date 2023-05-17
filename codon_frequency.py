import sys
from Bio import SeqIO
import pandas as pd

def count_letters(strings):
    letter_counts = []
    
    for string in strings:
        counts = {}
        for letter in string:
            if letter.isalpha():
                if letter in counts:
                    counts[letter] += 1
                else:
                    counts[letter] = 1
        letter_counts.append(counts)
    
    return letter_counts


def calculate_frequencies(letter_counts):
    frequencies = []
    
    for counts in letter_counts:
        total_letters = sum(counts.values())
        frequency = {letter: (count / total_letters)*100 for letter, count in counts.items()}
        frequencies.append(frequency)
    
    return frequencies


def count_letters_and_calculate_frequencies(sequences):
    letter_counts = count_letters(sequences)
    frequencies = calculate_frequencies(letter_counts)
    
    return frequencies


def convert_to_dataframe(frequencies, sequence_names):
    df = pd.DataFrame(frequencies)
    df.insert(0, "Sequence", sequence_names)
    df = df.fillna(0)  # Replace NaN values with 0
    return df


def visualize_letter_frequencies(fasta_file, output_file):
    sequences = []
    sequence_names = []
    
    # Read the FASTA file and extract sequences and names
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences.append(str(record.seq))
        sequence_names.append(record.id)
    
    frequencies = count_letters_and_calculate_frequencies(sequences)
    df = convert_to_dataframe(frequencies, sequence_names)
    
    # Write the dataframe to a CSV file
    df.to_csv(output_file, index=False)


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Please provide a FASTA file and an output file name as arguments.")
    else:
        fasta_file = sys.argv[1]
        output_file = sys.argv[2]
        visualize_letter_frequencies(fasta_file, output_file)
