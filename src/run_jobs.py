from orca_jobs import orca_job_sequence


if __name__ == "__main__":

    # Apparently, B3LYP is fine for geometry optimizations here.
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="6-311++G(d,p)",
        newgto='46 "LANL2DZ"',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )
    single_pt_arguments = dict(
        functional="CAM-B3LYP",
        basis_set="6-311++G(d,p)",
        newgto='46 "LANL2DZ"',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/dft_test',
        destination_path='data/batch_1_CAM_B3LYP',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )
