import pandas as pd
from Bio import SeqIO


def kidera_factors(input_file):
    """
    Compute Kidera factors for a given sequence of amino acids.

    Parameters:
         input_file : path of the fasta file

    """

    column_header = ['KF1_Helix/Bend_Preference', 'KF2_Side_Chain_Size', 'KF3_Extended_Structure_Preference',
              'KF4_Hydrophobicity', 'KF5_Double_Bend_Preference', 'KF6_Partial_Specific_Volume',
              'KF7_Flat_Extended_Preference', 'KF8_Occurence_In_Alpha_Region', 'KF9_pK-C',
              'KF10_Surrounding_Hydrophobicity']
    kidera_data = pd.DataFrame.from_records(list(map(lambda x: list(map(float, x.split(','))),
                                                     "-1.56,-1.67,-0.97,-0.27,-0.93,-0.78,-0.2,-0.08,0.21,-0.48;0.22,1.27,1.37,1.87,-1.7,0.46,0.92,-0.39,0.23,0.93;1.14,-0.07,-0.12,0.81,0.18,0.37,-0.09,1.23,1.1,-1.73;0.58,-0.22,-1.58,0.81,-0.92,0.15,-1.52,0.47,0.76,0.7;0.12,-0.89,0.45,-1.05,-0.71,2.41,1.52,-0.69,1.13,1.1;-0.47,0.24,0.07,1.1,1.1,0.59,0.84,-0.71,-0.03,-2.33;-1.45,0.19,-1.61,1.17,-1.31,0.4,0.04,0.38,-0.35,-0.12;1.46,-1.96,-0.23,-0.16,0.1,-0.11,1.32,2.36,-1.66,0.46;-0.41,0.52,-0.28,0.28,1.61,1.01,-1.85,0.47,1.13,1.63;-0.73,-0.16,1.79,-0.77,-0.54,0.03,-0.83,0.51,0.66,-1.78;-1.04,0,-0.24,-1.1,-0.55,-2.05,0.96,-0.76,0.45,0.93;-0.34,0.82,-0.23,1.7,1.54,-1.62,1.15,-0.08,-0.48,0.6;-1.4,0.18,-0.42,-0.73,2,1.52,0.26,0.11,-1.27,0.27;-0.21,0.98,-0.36,-1.43,0.22,-0.81,0.67,1.1,1.71,-0.44;2.06,-0.33,-1.15,-0.75,0.88,-0.45,0.3,-2.3,0.74,-0.28;0.81,-1.08,0.16,0.42,-0.21,-0.43,-1.89,-1.15,-0.97,-0.23;0.26,-0.7,1.21,0.63,-0.1,0.21,0.24,-1.15,-0.56,0.19;0.3,2.1,-0.72,-1.57,-1.16,0.57,-0.48,-0.4,-2.3,-0.6;1.38,1.48,0.8,-0.56,0,-0.68,-0.31,1.03,-0.05,0.53;-0.74,-0.71,2.04,-0.4,0.5,-0.81,-1.07,0.06,-0.46,0.65".split(
                                                         ";"))),
                                            index=["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F",
                                                   "P", "S", "T", "W", "Y", "V"], columns=column_header)

    records = SeqIO.parse(input_file, 'fasta')

    seq = []
    ids = []
    for record in records:
        seq.append(str(record.seq))
        ids.append(record.id)


    # Compute Kidera factors for each sequence
    k_factors = []
    for sequence in seq:
        factor_values = {}
        for factor in kidera_data.columns:
            # Compute the average of the factor values for the sequence
            factor_values[factor] = kidera_data.loc[list(sequence), factor].mean()
        k_factors.append(factor_values)

    # Convert the results to a DataFrame
    result_df = pd.DataFrame(k_factors)
    result_df
    result_df.insert(0, 'Protein_ID', ids)

    result_df.to_csv('kiderafactors_results.csv', index=False, header=True)

    


