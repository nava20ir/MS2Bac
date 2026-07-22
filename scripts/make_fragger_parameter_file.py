PARAM_TEMPLATE = """
database_name = /fragger_projects/{fasta_file}
num_threads = 1			# 0=poll CPU to set num threads; else specify num threads directly (max 64)

precursor_mass_lower = -10
precursor_mass_upper = 10
precursor_mass_units = 1			# 0=Daltons, 1=ppm, 3=DIA, 2=DIA_MS1
precursor_true_tolerance = 20
precursor_true_units = 1			# 0=Daltons, 1=ppm
fragment_mass_tolerance = 20
fragment_mass_units = 1			# 0=Daltons, 1=ppm
calibrate_mass = 0			# 0=Off, 1=On, 2=On and find optimal parameters
decoy_prefix = rev_

deisotope = 1
isotope_error = 0/1			# 0=off, -1/0/1/2/3 (standard C13 error)
mass_offsets = 0			# allow for additional precursor mass window shifts. Multiplexed with isotope_error. mass_offsets = 0/79.966 can be used as a restricted ‘open’ search that looks for unmodified and phosphorylated peptides (on any residue)
precursor_mass_mode = selected

remove_precursor_peak = 0			# 0 = not remove, 1 = only remove the peak with the precursor charge, 2 = remove all peaks with all charge states. Default: 0
remove_precursor_range = -1.50,1.50			# Unit: Da. Default: -1.5,1.5
intensity_transform = 0			# 0 = none, 1 = sqrt root. Default: 0

write_calibrated_mgf = 0			# Write calibrated MS2 scan to a MGF file (0 for No, 1 for Yes).
mass_diff_to_variable_mod = 0			# Put mass diff as a variable modification. 0 for no; 1 for yes and change the original mass diff and the calculated mass accordingly; 2 for yes but not change the original mass diff and the calculated mass.

localize_delta_mass = 0
delta_mass_exclude_ranges = (-1.5,3.5)
fragment_ion_series = b,y
ion_series_definitions = 

labile_search_mode = off			# options: "off", "nglycan", "labile" (default: off)
deltamass_allowed_residues = ST			# aminoacids that are allowed to be modified in Glyco mode. E.g. "ST"
diagnostic_intensity_filter = 0			# possible values are 0 <= x <= 1
Y_type_masses = 0/203.07937/406.15874/568.21156/730.26438/892.3172/349.137279
diagnostic_fragments = 204.086646/186.076086/168.065526/366.139466/144.0656/138.055/126.055/163.060096/512.197375/292.1026925/274.0921325/657.2349/243.026426/405.079246/485.045576/308.09761

search_enzyme_name = NoCut
search_enzyme_cutafter = Z
num_enzyme_termini = 2			# 2 for enzymatic, 1 for semi-enzymatic, 0 for nonspecific digestion
allowed_missed_cleavage = 0			# maximum value is 5
clip_nTerm_M = 0

#maximum of 7 mods - amino acid codes, * for any amino acid, [ and ] specifies protein termini, n and c specifies peptide termini
# variable_mod_01 = 15.99490 M 1
# variable_mod_02 = 42.01060 [^ 1
# variable_mod_03 = 79.96633 STY 3
# variable_mod_04 = -17.02650 nQnC 1
# variable_mod_05 = -18.01060 nE 1
# variable_mod_06 = 0.00000 site_06 1
# variable_mod_07 = 0.00000 site_07 1
# variable_mod_08 = 0.00000 site_08 1
# variable_mod_09 = 0.00000 site_09 1
# variable_mod_10 = 0.00000 site_10 1
# variable_mod_11 = 0.00000 site_11 1
# variable_mod_12 = 0.00000 site_12 1
# variable_mod_13 = 0.00000 site_13 1
# variable_mod_14 = 0.00000 site_14 1
# variable_mod_15 = 0.00000 site_15 1
# variable_mod_16 = 0.00000 site_16 1

allow_multiple_variable_mods_on_residue = 0			# static mods are not considered
max_variable_mods_per_mod = 3
max_variable_mods_per_peptide = 0			# maximum 5
max_variable_mods_combinations = 5000			# maximum 65534, limits number of modified peptides generated from sequence

output_file_extension = pepXML
output_format = tsv_pepXML
output_report_topN = 1
output_max_expect = 50
report_alternative_proteins = 1			# 0=no, 1=yes

precursor_charge = 1 4			# precursor charge range to analyze; does not override any existing charge; 0 as 1st entry ignores parameter
override_charge = 0			# 0=no, 1=yes to override existing precursor charge states with precursor_charge parameter

digest_min_length = 7
digest_max_length = 30
digest_mass_range = 700.0 5000.0			# MH+ peptide mass range to analyze
max_fragment_charge = 2			# set maximum fragment charge state to analyze (allowed max 5)

#open search parameters
track_zero_topN = 0			# in addition to topN results, keep track of top results in zero bin
zero_bin_accept_expect = 0.00			# boost top zero bin entry to top if it has expect under 0.01 - set to 0 to disable
zero_bin_mult_expect = 1.00			# disabled if above passes - multiply expect of zero bin for ordering purposes (does not affect reported expect)
add_topN_complementary = 0

# spectral processing

minimum_peaks = 15			# required minimum number of peaks in spectrum to search (default 10)
use_topN_peaks = 500
min_fragments_modelling = 2
min_matched_fragments = 4
minimum_ratio = 0.01			# filter peaks below this fraction of strongest peak
clear_mz_range = 0.0 0.0			# for iTRAQ/TMT type data; will clear out all peaks in the specified m/z range

# additional modifications

add_Cterm_peptide = 0.000000
add_Nterm_peptide = 0.000000
add_Cterm_protein = 0.000000
add_Nterm_protein = 0.000000

add_G_glycine = 0.000000
add_A_alanine = 0.000000
add_S_serine = 0.000000
add_P_proline = 0.000000
add_V_valine = 0.000000
add_T_threonine = 0.000000
add_C_cysteine = 57.021460
add_L_leucine = 0.000000
add_I_isoleucine = 0.000000
add_N_asparagine = 0.000000
add_D_aspartic_acid = 0.000000
add_Q_glutamine = 0.000000
add_K_lysine = 0.000000
add_E_glutamic_acid = 0.000000
add_M_methionine = 0.000000
add_H_histidine = 0.000000
add_F_phenylalanine = 0.000000
add_R_arginine = 0.000000
add_Y_tyrosine = 0.000000
add_W_tryptophan = 0.000000
add_B_user_amino_acid = 0.000000
add_J_user_amino_acid = 0.000000
add_O_user_amino_acid = 0.000000
add_U_user_amino_acid = 0.000000
add_X_user_amino_acid = 0.000000
add_Z_user_amino_acid = 0.000000
"""


from pathlib import Path


def write_param_file(fasta_file, output_file=None):
    

    if output_file is None:
        output_file = fasta_file.with_suffix(".params")

    with open(output_file, "w") as f:
        f.write(
            PARAM_TEMPLATE.format(
                fasta_file=fasta_file
            )
        )

    return output_file

if __name__ == '__main__':
    import sys
    write_param_file(sys.argv[1], output_file=sys.argv[2])
