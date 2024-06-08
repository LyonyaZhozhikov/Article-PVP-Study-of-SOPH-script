import os
from datetime import datetime
import Bio.PDB

while True:
    protein_name = input("Enter protein name (or 'q' to quit): ")
    if protein_name == 'q':
        break
    amino_acids = input("Enter whole amino acid sequence: ")
    while True:
        # timestamp = datetime.now().strftime("%H-%M-%S")
        folder_date = datetime.now().strftime("%m-%d-%Y")
        folder_name = f"{protein_name}_{folder_date}"
        os.makedirs(f'proteins_esm/{folder_name}', exist_ok=True)

        # separation
        position_input = input("Enter position and variant like this - '254D,1914H and so on' (or 'n' for new protein): ")
        if position_input == 'n':
            break
        aa_numbers = []
        aa_change = []

        for item in position_input.split(','):
            aa_numbers.append(int(''.join(filter(str.isdigit, item))))
            aa_change.append(''.join(filter(str.isalpha, item)))
        # Prompt the user to enter the step size and length
        step_size = int(input("Enter the step size for RMSD calculation: "))
        step_length = int(input("Enter the step length for RMSD calculation: "))

        # Calculate the new start and end positions based on the step size and length
        rmsd_start_cut = list(range(0, len(amino_acids), step_size))
        rmsd_end_cut = [start + step_length for start in rmsd_start_cut]
        for i in range(len(aa_numbers)):
            # RMSD part
            # cycle
            position = int(aa_numbers[i])
            substitution = str(aa_change[i])
            for j in range(len(rmsd_start_cut)):
                start_id = rmsd_start_cut[j]
                end_id = rmsd_end_cut[j]
                atoms_to_be_aligned = range(start_id, end_id + 1)

                # Start the parser
                pdb_parser = Bio.PDB.PDBParser(QUIET=True)

                # Get the structures
                ref_structure = pdb_parser.get_structure("reference", f"proteins_esm/{folder_name}/{protein_name}_{position}.pdb")
                sample_structure = pdb_parser.get_structure("sample", f"proteins_esm/{folder_name}/{protein_name}_{position}_{substitution}_mutated.pdb")

                # Use the first model in the pdb-files for alignment
                # Change the number 0 if you want to align to another structure
                ref_model = ref_structure[0]
                sample_model = sample_structure[0]

                # Make a list of the atoms (in the structures) you wish to align.
                # In this case we use CA atoms whose index is in the specified range
                ref_atoms = []
                sample_atoms = []

                # Iterate of all chains in the model in order to find all residues
                for ref_chain in ref_model:
                    # Iterate of all residues in each model in order to find proper atoms
                    for ref_res in ref_chain:
                        # Check if residue number ( .get_id() ) is in the list
                        if ref_res.get_id()[1] in atoms_to_be_aligned:
                            # Append CA atom to list
                            ref_atoms.append(ref_res['CA'])

                # Do the same for the sample structure
                for sample_chain in sample_model:
                    for sample_res in sample_chain:
                        if sample_res.get_id()[1] in atoms_to_be_aligned:
                            sample_atoms.append(sample_res['CA'])

                # Now we initiate the superimposer:
                super_imposer = Bio.PDB.Superimposer()
                super_imposer.set_atoms(ref_atoms, sample_atoms)
                super_imposer.apply(sample_model.get_atoms())

                # Print RMSD:
                rmsd = super_imposer.rms

                def save_rmsd_to_file(rmsd_s, protein_name_s, position_s, start_id_s, end_id_s):
                    start_position = position_s + start_id_s - 200
                    end_position = position_s + end_id_s - 200
                    with open(f"RMSD/{protein_name_s}_rmsd_log.txt", "a") as f:
                        f.write(f"RMSD at position {start_position} {end_position} - {rmsd_s}\n")

                save_rmsd_to_file(rmsd, protein_name, position, start_id, end_id)
        #
        with open(f"RMSD/{protein_name}_rmsd_log_1.txt", "a") as f:
            f.write(f"\n")
            # Save the aligned version of 1UBQ.pdb
            # io = Bio.PDB.PDBIO()
            # io.set_structure(sample_structure)
            # io.save("NBAS_07-17-2023/NBAS_1914_aligned.pdb")
