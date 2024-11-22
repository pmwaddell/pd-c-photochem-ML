from orca_jobs import orca_job_sequence


if __name__ == "__main__":
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVP",
        RI="RIJCOSX",
        dispersion_correction="D3BJ"
    )
    single_pt_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVP",
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/dft_test',
        destination_path='data/batch_1_B3LYP',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )
