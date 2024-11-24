import logging
logger = logging.getLogger(__name__)

from orca_jobs import orca_job_sequence


def configure_logging(log_filename: str="logs/log.log") -> None:
    """Configure the logger for our purposes."""
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s %(message)s',
        datefmt='%H:%M:%S',
        handlers=[logging.FileHandler(f"logs/{log_filename}.log", mode='w'),
                  stream_handler]
    )


if __name__ == "__main__":
    configure_logging("lot_benchmarks.log")

    # TODO: seems to be giving us problems with Pd, might have to work out newgto again?
    # TODO skip for now but revisit this...
    # # 6311
    # geom_opt_arguments = dict(
    #     functional="B3LYP",
    #     basis_set="6-311++G(d,p)",
    #     newgto='',
    #     RI="RIJCOSX",
    #     dispersion_correction="D3BJ",
    #     solvent="CPCM(Chloroform)"
    # )
    #
    # single_pt_arguments = dict(
    #     functional="CAM-B3LYP",
    #     basis_set="6-311++G(d,p)",
    #     newgto='',
    #     RI="RIJCOSX",
    #     dispersion_correction="D3BJ",
    #     solvent="CPCM(Chloroform)",
    #     NMR=True,
    #     freq=True
    # )
    #
    # orca_job_sequence(
    #     path_to_conf_search_xyz_files='data/LOT_benchmarks/LOT_benchmark_inputs',
    #     destination_path='data/LOT_benchmarks/basis_sets/6311',
    #     geom_opt_arguments=geom_opt_arguments,
    #     single_pt_arguments=single_pt_arguments
    # )

    # def2-TZVP
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )

    single_pt_arguments = dict(
        functional="CAM-B3LYP",
        basis_set="def2-TZVP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/LOT_benchmarks/LOT_benchmark_inputs',
        destination_path='data/LOT_benchmarks/basis_sets/def2-TZVP',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )

    # def2-TZVPP
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVPP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )

    single_pt_arguments = dict(
        functional="CAM-B3LYP",
        basis_set="def2-TZVPP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/LOT_benchmarks/LOT_benchmark_inputs',
        destination_path='data/LOT_benchmarks/basis_sets/def2-TZVPP',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )

    # def2-TZVPD
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVPD",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )

    single_pt_arguments = dict(
        functional="CAM-B3LYP",
        basis_set="def2-TZVPD",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/LOT_benchmarks/LOT_benchmark_inputs',
        destination_path='data/LOT_benchmarks/basis_sets/def2-TZVPD',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )

    # pcseg-2
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="pcseg-2",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )

    single_pt_arguments = dict(
        functional="CAM-B3LYP",
        basis_set="pcseg-2",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/LOT_benchmarks/LOT_benchmark_inputs',
        destination_path='data/LOT_benchmarks/basis_sets/pcseg-2',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )

    # aug-pcseg-2
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="aug-pcseg-2",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)"
    )

    single_pt_arguments = dict(
        functional="CAM-B3LYP",
        basis_set="aug-pcseg-2",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="CPCM(Chloroform)",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/LOT_benchmarks/LOT_benchmark_inputs',
        destination_path='data/LOT_benchmarks/basis_sets/aug-pcseg-2',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )

    # def2-TZVP NO SOLVENT
    geom_opt_arguments = dict(
        functional="B3LYP",
        basis_set="def2-TZVP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent=""
    )

    single_pt_arguments = dict(
        functional="CAM-B3LYP",
        basis_set="def2-TZVP",
        newgto='',
        RI="RIJCOSX",
        dispersion_correction="D3BJ",
        solvent="",
        NMR=True,
        freq=True
    )

    orca_job_sequence(
        path_to_conf_search_xyz_files='data/LOT_benchmarks/LOT_benchmark_inputs',
        destination_path='data/LOT_benchmarks/basis_sets/def2-TZVP_no_solvent',
        geom_opt_arguments=geom_opt_arguments,
        single_pt_arguments=single_pt_arguments
    )